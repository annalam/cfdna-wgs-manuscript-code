use parse_args;
use std::str;
use std::io::{Write, BufReader, BufRead, stderr};
use std::process::exit;
use std::collections::VecDeque;
use std::path::Path;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::{Cigar, Record};

use genome_reader::RefGenomeReader;
use reader_of_bams::ReaderOfBams;

const MAX_READ_LEN: usize = 5000;
const WINDOW_LEN: usize = 50000;
const USAGE: &str = "
Usage:
  mutato call [options] <genome.fa> <bam_files>...
  mutato call2 [options] <genome.fa> <bam_files>...

Options:
  --region=REGION     Genomic region to include in analysis [default: all]
  --alt-reads=N       Minimum alt allele reads to report [default: 0]
  --alt-frac=N        Minimum alt allele fraction to report [default: 0.0]
  --min-mapq=N        Ignore reads with mapping quality < N [default: 0]
  --min-baseq=N       Ignore bases with a base quality < N [default: 0]
  --count-duplicates  Count reads that are marked as duplicates
  --max-frag-len=N    Ignore DNA fragments longer than N bases [default: Inf]
  --verbose           More output during processing 
";


// Information from the reads are collected into "pileups" for each
// genomic coordinate. The pileups are processed for windows of the
// genome.

#[derive(Copy, Clone)]
struct Substitution {
	reads: u32,
	sidedness: f32, // Average position from nearest read end
	mapq: f32
}

impl Substitution {
	pub fn new() -> Substitution {
		Substitution { reads: 0, sidedness: 0.0, mapq: 0.0 }
	}
}

#[derive(Clone)]
struct Indel {
	sequence: Box<str>,    // Prefixed with '+' for insertion, '-' for deletion
	reads: u32,
	sidedness: f32,
	mapq: f32
}


// TODO: We could merge Substitution and Indel into one Mutation type with
// a 64-bit field that specifies the alternate allele. Then we could have
// Pileup simply be a SmallVec with room for 1 or 2 Mutations objects.
#[derive(Clone)]
struct Pileup {
	total: u32,
	acgt: [Substitution; 4],
	indels: Vec<Indel>
}


fn reset_pileups( pileups : &mut Vec<VecDeque<Pileup>>, n_files : usize, buffer_len : usize )
{
		pileups.clear();
		// Generate new pileup vectors
		for b in 0..n_files {
			pileups.push( VecDeque::with_capacity( buffer_len));
			fill_pileup( &mut pileups[ b], buffer_len);
		}
}

fn fill_pileup( pileup : &mut VecDeque<Pileup>, buffer_len : usize ){

	for _k in 0..(buffer_len - pileup.len()) {

		pileup.push_back( Pileup { total: 0, acgt: [
						Substitution::new(), Substitution::new(),
						Substitution::new(), Substitution::new()
						], indels: Vec::new() });
	}
}

pub fn main() {

	let args = parse_args( USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let bam_paths = args.get_vec("<bam_files>");
	let alt_reads: u32 = args.get_str("--alt-reads").parse().unwrap_or_else(
		|_| error!("--alt-reads must be a positive integer"));
	let alt_frac: f32 = args.get_str("--alt-frac").parse().unwrap_or_else(
		|_| error!("--alt-frac must be a fraction between 0 and 1."));
	let min_mapq: u8 = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be a non-negative integer."));
	let min_baseq: u8 = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a non-negative integer."));
	let count_duplicates = args.get_bool("--count-duplicates");
	let max_frag_len: usize = if args.get_str("--max-frag-len") == "Inf" {
		0
	} else {
		args.get_str("--max-frag-len").parse().unwrap()
	};
	let genome_region = args.get_str("--region");
	let verbose = args.get_bool("--verbose");

	// Read the reference genome FASTA index into memory.
	let mut ref_genome_reader = RefGenomeReader::new(genome_path);	

	let n_bams = bam_paths.len();
	let mut pileups: Vec<VecDeque<Pileup>> = Vec::with_capacity(n_bams);

	let mut bam_reader = ReaderOfBams::new(n_bams);

	// Set up bam file reader that makes sure
	// all input bam files are congruent (consistent names and lengths)
	// ReaderOfBams uses indexed reader if index files are available
	// If index files are not available it iterates through all
	// records
	for b in 0..n_bams { bam_reader.add_file( &bam_paths[ b]); }

	//if bam_reader.chr_names.contains( String::from( &genome_region)) { region_is_valid = true; }
	let region_is_valid = bam_reader.index_for_name(genome_region).is_ok();
	let mut next_chr = 0i32;

	if !region_is_valid && genome_region != "" && genome_region != "all" {
		error!("Invalid genomic region '{}' specified.", genome_region);
	} else if genome_region != "" && genome_region != "all" {
		next_chr = bam_reader.index_for_name(genome_region)
			.unwrap_or_else(|_| error!("Cannot find index for region '{}'", genome_region)) as i32;
	}

	// DEBUG - Start from chr index
	//next_chr = 454;

	// Use path from first bam file
	//let path = Path::new( bam_paths[ 0]).parent().unwrap();

	// Output VCF file information
	//println!("##fileformat=nonstandard_VCF");
	//println!("##source=Mutato call2_V0.9");
	//println!("##FORMAT=<ID=AL,Number=1,Type=Integer,Description=\"Number of alternate reads\">");
	//println!("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
	//println!("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality\">");
	//println!("##FORMAT=<ID=BQ,Number=1,Type=Integer,Description=\"Base quality\">");
	//println!("##FORMAT=<ID=SD,Number=1,Type=Integer,Description=\"Read sidedness\">");
	//println!("##BAM_BASE_URL={}", path.to_str().unwrap());

	// Output header row.
	print!("CHROM\tPOSITION\tREF\tALT\tNOTES");
	for bam_path in &bam_paths {
		// Remove bam-suffix and possible filepath
		let path = Path::new( bam_path);
		print!("\t{}", path.file_stem().unwrap().to_str().unwrap());
	}
	println!();

	let buffer_len = WINDOW_LEN + MAX_READ_LEN;

	let mut skip_from_chr = -1i32; 	

	let mut loaded_ref_chr = -1;
	let mut chr_seq = Vec::new();

	eprintln!( "INFO: Running mutato call2...");

	// Loop through all regions
	'chr_loop: for chr_index in 0..bam_reader.chr_names.len() {
		
		let chr_index_i32 = chr_index as i32;
		// Skipping chromosomes
		if next_chr >= 0 && next_chr > chr_index_i32 { 
			if skip_from_chr == -1 { skip_from_chr = chr_index_i32; }
			continue; 
		}

		// Print skipping info when skipped
		if skip_from_chr >= 0 {
			if verbose { eprintln!( "INFO: Skipped from region index: {} -> {}", skip_from_chr, chr_index); }
			skip_from_chr = -1;
		}
		next_chr = -1;
		
		// Skipping a region will cause the records to be read and discarded
		// Genomic regions (chromosomes) are skipped if the --region flag is used
		let skip_region = if genome_region == "" || genome_region == "all" { false } 
						  else if genome_region == bam_reader.chr_names[ chr_index]  { false }
						  else { true };

		if skip_region && verbose {
			eprintln!( "INFO: Skipping {}...", &bam_reader.chr_names[ chr_index]);
		}

	    // If we are skipping a region after the region has already been
	    // found we are done!	   
		if skip_region { 
			eprintln!( "INFO: Region {} processed. Exitting...", &genome_region);
			break 'chr_loop; 
		}		
		
		// Reset the pileup vectors, because we are in a new chromosome.
		reset_pileups( &mut pileups, n_bams, buffer_len);	

		// These two variables keep track of the current pileup window
		// Inside the current reference chromosome (region)
		let mut window_start;
		let mut window_end: usize = 0;	
		//let mut window_next = 0usize;	

		let mut last_window = false; // Switching region (chromosome) after this window		
		let region_len = bam_reader.chr_lens[ chr_index]; 
		let mut read_reach: usize = 0;
		
		if verbose { eprintln!("INFO: CHR INDEX: {}/{} ({})", chr_index, bam_reader.chr_names.len()-1, bam_reader.chr_names[ chr_index]); }

		// Loop while reads inside the read window
		// TODO: Use for-loop with .step_by() here when stable
		'window_loop: loop {  	

			// Move window to next position
			window_start = window_end;
			window_end = window_start + WINDOW_LEN;

			if window_start > bam_reader.chr_lens[ chr_index] {
				error!("Window set beyond region '{}' length {} at {}.",bam_reader.chr_names[ chr_index], region_len, window_start);
			} else if window_end >= bam_reader.chr_lens[ chr_index] { last_window = true;
			} // if window reaches beyond end of chr, it is the last window

			// DEBUG window tracking
			if verbose {
				eprintln!( "WINDOW: {} {}-{} ({:.2}%) {}",
					bam_reader.chr_names[ chr_index],
					window_start, window_end,
					(window_start as f32 / (region_len+1)as f32) * 100.0,
					if last_window {&"LAST"} else {&""}  );
			}
			
			// Go through the window in every input BAM file.
			'bam_loop: for bam_index in 0..n_bams {
			
				// Append enough empty Pileup structures into the (drained) buffer,
				// so we do not try to write past its boundaries.
				fill_pileup( &mut pileups[ bam_index], buffer_len);

				'reads_loop: loop {

					// read_next_record returns false if there are no more reads
					// for the specified region and window
					// Use window end +1 to be able to assign indels at the window 
					// border to the previous genomic coordinate
					if bam_reader.read_next_record( bam_index, chr_index_i32, window_start, window_end+1 ) == false {
						// No more records for the requested window
						// Move on to next bam file
						break 'reads_loop; 
					}		

					//New record (read) was loaded
					let ref read = bam_reader.record.as_ref()
						.unwrap_or_else(|| error!("Record not set."));

					let pos = read.pos() as usize;    // 0-based position	
					let tid = read.tid();  // chr index, -1 if none

					// Sanity checks
					if pos < window_start && tid == chr_index_i32 { 
						eprintln!("Window set to: {}:{} {}-{}",bam_reader.chr_names[ chr_index], chr_index, window_start, window_end);
						error!("Reads in BAM file '{}' are not sorted by position. ({}:{}-{})", bam_paths[ bam_index], tid, pos, window_start); 
					} else if tid != chr_index_i32 {
						error!("Reader returned record from incorrect region. chr:{} vs. read:{}", tid, chr_index); 
					}

					// Check if the read needs to be processed
					if read.is_unmapped() { continue; }
					if read.is_secondary() { continue; }
					if read.is_supplementary() { continue; }
					if read.is_quality_check_failed() { continue; }

					// User can choose to count reads that
					// have been flagged as PCR/optical duplicates.
					if read.is_duplicate() && !count_duplicates { continue; }

					if read.mapq() < min_mapq { continue; }

					// If the user only wants to count concordant reads, check
					// that the mates have aligned to the same chromosome and
					// are not too far apart.
					if max_frag_len > 0 {
						if read.tid() != read.mtid() { continue; }
						if read.insert_size().abs() > max_frag_len as i64 {
							continue;
						}
					}

					// Load required reference sequence for the read,
					// if not already loaded
					if loaded_ref_chr != tid { 
						chr_seq = ref_genome_reader.load_chromosome_seq( &bam_reader.chr_names[ chr_index]); 
						loaded_ref_chr = tid;
					}														
								
					// Process read CIGAR and add its information to pileups
					let seq = read.seq();
					let cigar = read.cigar();
					let qual = read.qual();
					let read_length = seq.len();
					if read_length > MAX_READ_LEN {
						eprintln!("WARNING: Read {} is {} bases long, exceeding the maximum length threshold of {}.", str::from_utf8(read.qname()).unwrap(), seq.len(), MAX_READ_LEN);
						continue;
					}

					read_reach = ::std::cmp::max(read_reach, pos + read_length);

					// These variables keep track of our current position
					// within the pileup track and the read sequence, as we
					// read through the CIGAR string.
					let mut pileup_idx = pos - window_start;
					let mut seq_idx = 0;
					let mut prior_n = false;

					// Handle all CIGAR elements in order
					'cigar_loop: for s in cigar.iter() { match *s {

						Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
							for _ in 0..len {
								let acgt_idx = encoded_base_to_acgt_index( seq.encoded_base(seq_idx));
								if acgt_idx >= 4 {
									prior_n = true;  // Ambiguous base
								} else if qual[seq_idx] < min_baseq {
									prior_n = true;  // Low quality base, ignore
								} else {
									let mut pileup = &mut pileups[ bam_index][ pileup_idx];
									let mut acgt = &mut pileup.acgt[ acgt_idx];

									prior_n = false;
									pileup.total += 1;
									
									acgt.reads += 1;
									let n_reads = acgt.reads as f32;
									let n_reads_before = n_reads - 1.0;
									
									// Distance from nearest read end
									let end_distance = ::std::cmp::min( seq_idx, (seq.len()-1)-seq_idx ) as f32;

									// Calculate cumulative moving averages
									acgt.sidedness = (acgt.sidedness * n_reads_before + end_distance) / n_reads;
									acgt.mapq = (acgt.mapq * n_reads_before + read.mapq() as f32) / n_reads;
								}
								pileup_idx += 1;
								seq_idx += 1;
							}
						},
						Cigar::Ins(len) => {

							if pileup_idx == 0 || seq_idx == 0 { 
								eprintln!( "WARNING: Insertion at index 0 cannot be assigned to correct genome position." );
								seq_idx += len as usize;
								continue;
							}
					
							let mut allele = String::with_capacity(len as usize + 2);
							allele.push('+'); //Mark insertions with a plus
							allele.push( (chr_seq[ window_start + pileup_idx-1] as char).to_ascii_uppercase());

							for l in 0..len {
								allele.push( encoded_base_to_char( seq.encoded_base( seq_idx)));
								//base_qual_avg = (base_qual_avg * (l as f32) + (qual[ seq_idx] as f32)) / (l+1) as f32;
								seq_idx += 1 as usize;
							}

							// Discard insertions with any Ns
							if allele.contains("N") { continue; } 

							// Distance from nearest read end
							let end_distance = ::std::cmp::min( seq_idx-(len as usize), (seq.len()-1)-seq_idx ) as f32;
							let pileup = &mut pileups[ bam_index][ pileup_idx-1];

							count_indel( pileup, &allele, end_distance, read.mapq() as f32);
							// Indels are assigned to the coordinate of the previous MATCHing base
							// If the previous MATCHing base was an N, it has not been added to 
							// pile.total and needs to be added here
							if prior_n { pileup.total += 1; }
						},
						Cigar::Del(len) => {

							if pileup_idx == 0 { 
								eprintln!( "WARNING: Deletion at index 0 cannot be assigned to correct genome position." );
								pileup_idx += len as usize;
								continue;
							}
							//let mut pileup = &mut pileups[ bam_index].get_mut( pileup_idx).unwrap();
							let mut pileup = &mut pileups[ bam_index][ pileup_idx-1];							
							let mut allele = String::with_capacity( len as usize + 2);
							allele.push('-');
							allele.push( (chr_seq[ window_start + pileup_idx-1] as char).to_ascii_uppercase());

							for _ in 0..len {
								allele.push( (chr_seq[ window_start + pileup_idx] as char).to_ascii_uppercase());
								pileup_idx += 1;
							}

							// No N in ref genome
							// if allele.contains("N") { continue; }

							// Distance from nearest read end
							let end_distance = ::std::cmp::min( seq_idx, (seq.len()-1)-seq_idx ) as f32;								
							count_indel( &mut pileup, &allele, end_distance, read.mapq() as f32);
							// Indels are assigned to the coordinate of the previous MATCHing base
							// If the previous MATCHing base was an N, it has not been added to 
							// pile.total and needs to be added here								
							if prior_n { pileup.total += 1; }
						},
						Cigar::SoftClip(_len) | Cigar::HardClip(_len) => {
							// Soft and hard clips must be the last item
							// in a CIGAR string. So if we see them,
							// we are done with the read.
							break 'cigar_loop;
						},
						Cigar::RefSkip(_len) =>
							error!("Unexpected CIGAR type: N"),
						Cigar::Pad(_len) =>
							error!("Unexpected CIGAR type: P")
					}} //for each cigar
				} // for each read
			} //for each bam
		
			// Process pileups if there are any reads within the window.	
			// If no reads begin or extend into the current window,
			// pileups are empty
			if read_reach > window_start { 

				assert!(loaded_ref_chr >= 0);
				//eprintln!( "processing window: {}, {}-{} (reglen:{})", chr_index, window_start, window_end, region_len);

				// Pileups have been filled for the window + buffer
				// Output variations that are above specified tresholds
				process_pileups( &pileups, chr_index, window_start, 0, WINDOW_LEN, &chr_seq, alt_frac, alt_reads, &bam_reader );
				
				// Pileups will be reset after the last window
				// and do not need to be drained
				if !last_window {
					// eprintln!("INFO: Draining window pileups");
					for b in 0..n_bams {
						// Drain (delete) one window's worth of positions from the buffer
						pileups[ b].drain(0..WINDOW_LEN);
					}
				}
			}		

			// After last window in the region, 
			// move on to next region
			if last_window { break 'window_loop; }

		} // Window inside chromosome
	} // For each chromosome
}

// Print pileups for genomic positions where at least one sample
// has sufficient evidence of variation
// NOTE: start_ind and end_ind are relative to window_begin
fn process_pileups( pileups : &Vec<VecDeque<Pileup>>, chr_index : usize, window_begin : usize, start_ind : usize, end_ind : usize, ref_seq : &Vec<u8>,
                    alt_frac : f32, alt_reads : u32, bam_reader : &ReaderOfBams )  {

	let n_bams = pileups.len();
	let region_len = bam_reader.chr_lens[ chr_index];
	if start_ind >= end_ind { return; } 

	for k in start_ind..end_ind {
	
		//t_build_indels.start();				
		let seq_idx = window_begin + k;
	
		// Stop if we have reached the end of the region (chromosome).
		if seq_idx == region_len { break; }
	
		//let curpos = window_start + k;
		//if curpos % 1000000 == 0 { eprintln!( "POS: {} {}", &chr_names[ chr_index], curpos); }
		let ref_base = (ref_seq[ seq_idx] as char).to_ascii_uppercase();
	
		// Build a list of indels that need to be reported
		let mut report_indels: Vec<Box<str>> = Vec::new();
		for bam_index in 0..n_bams {
			
			let ref pileup = &pileups[ bam_index][ k];
	
			for i in 0..pileup.indels.len() {  
	
				let mut indel = &pileup.indels[ i];
				if (indel.reads as f32) / (pileup.total as f32) < alt_frac || indel.reads < alt_reads {
					continue; 
				}
	
				if report_indels.contains( &indel.sequence) == false {
					report_indels.push( indel.sequence.clone());
				}
			}
		}
	
		//t_build_indels.end();
		//t_print_indels.start();
	
		// Print the indels. Note that indels must be printed before
		// base substitutions, since their position will be pos-1, and
		// we want position-sorted output rows.
		for indel_seq in report_indels {
			//let prev_ref = &indel.sequence.chars().nth(0).unwrap(); 
			//let prev_ref = 'N';
			//let pref_ref = seq.encoded_base( seq_idx);
			let sign = indel_seq.chars().next().unwrap();
	
			if sign == '+' {
				print!("{}\t{}\t{}\t{}\t", bam_reader.chr_names[ chr_index], seq_idx+1, ref_base, &indel_seq[1..]);
			} else if sign == '-' {
				print!("{}\t{}\t{}\t{}\t", bam_reader.chr_names[ chr_index], seq_idx+1, &indel_seq[1..], ref_base);
			}
			else {
				error!("Unknown indel type.");
			}
	
			'outer: for j in 0..n_bams {
				let pile = pileups[ j].get( k).unwrap();
				for indel in &pile.indels {
					if indel.sequence == indel_seq {
						// Output indel stats
						print!("\t{}:{}:{:.0}::{:.0}", indel.reads, pile.total, indel.mapq, indel.sidedness );
						//print!("\t{}:{}", indel.reads, pile.total);

						// TODO: Lower these thresholds because the mutation is likely to be true positive
						if indel.reads as f32 / pile.total as f32 >= alt_frac && indel.reads >= alt_reads {
							print!(":*");
						}
						continue 'outer;
					}
				}
				print!("\t0:{}:0:0:0", pile.total);
			}
			println!();
		}
	
		//t_print_indels.end();
		//t_print_bs.start(); //"Print base substitutions"
			
		// Print base substitutions that are above thresholds.
		const ACGT: [char; 4] = ['A', 'C', 'G', 'T'];
		for base in 0..4 {
			if ref_base == ACGT[ base] { continue; }
	
			let mut variation_found = false;
	
			// Any sample above thresholds?
			for bam_index in 0..n_bams {
				let ref pile = &pileups[ bam_index][ k];
				let reads = pile.acgt[ base].reads;
				if reads >= alt_reads && (reads as f32) / (pile.total as f32) >= alt_frac { 
					variation_found = true;
					break; 
				}
			}
	
			if !variation_found { continue; }
	
			print!("{}\t{}\t{}\t{}\t", bam_reader.chr_names[ chr_index], seq_idx+1, ref_base, ACGT[ base]);
			for bam_index in 0..n_bams {
				let ref pile = &pileups[ bam_index][ k];
				let ref acgt = &pile.acgt[ base];
				let reads = acgt.reads;						
				//if reads < alt_reads || (reads as f32) / (pileup.total as f32) < alt_frac { continue; }

				// Output mismatch statistics							
				print!("\t{}:{}:{:.0}::{:.0}", reads, pile.total, acgt.mapq, acgt.sidedness);
				//print!("\t{}:{}", reads, pile.total);

				let cur_frac = (reads as f32) / (pile.total as f32);
				// TODO: Lower these thresholds because the mutation is likely to be true positive
				if cur_frac >= alt_frac && reads >= alt_reads {
					print!(":*");
					//printed.push(format!("{}", ACGT[ base]));
				}
			}
			println!();					
		}
		//t_print_bs.end();				
		
	} //processing pileup window
}


fn count_indel(pileup: &mut Pileup, seq: &str, end_distance: f32, mapq: f32) {
	for i in 0..pileup.indels.len() {
		let mut indel = &mut pileup.indels[ i];
		if &*indel.sequence == seq {
			let n_reads_before = indel.reads as f32;
			let n_reads = n_reads_before + 1.0;
			indel.reads += 1;			
			// Calculate cumulative averages
			// TODO: Could just count sums here, and do the division when
			// printing.
			indel.sidedness = (indel.sidedness * n_reads_before + end_distance) / n_reads;
			indel.mapq = (indel.mapq * n_reads_before + mapq) / n_reads;
			return;
		}
	}
	pileup.indels.push(Indel {
		sequence: seq.into(), 
		reads: 1, 
		sidedness: end_distance, 
		mapq: mapq
	});
}



fn encoded_base_to_acgt_index(encoded: u8) -> usize {
	match encoded {
		1 => 0,   // A
		2 => 1,   // C
		4 => 2,   // G
		8 => 3,   // T
		_ => 4    // N
	}
}

fn encoded_base_to_char(encoded: u8) -> char {
	match encoded { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', _ => 'N' }
}
