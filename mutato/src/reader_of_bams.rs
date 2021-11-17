use std::{str, fs};

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::{Cigar, Record};

// Reader of bams reads records from BAM files and works
// both with (Indexed mode) and without bam.bai index files
// (iterative mode). "region" and "chromosome" used as synonyms

// NOTE: Can be only used with bam files with positionally 
//       sorted records.
// NOTE: Processing unaligned reads with target id (tid) -1
//       is not supported

// Usage:
//  let mut bam_reader = ReaderOfBams::new( n_bams, ignore_bam_index_files);
//  for b in 0..n_bams { bam_reader.add_file( &bam_paths[ b]); }
// Loop until read_next_record returns false
//  if bam_reader.read_next_record( bam_index, chr_index_i32, window_start, window_end+1 ) == false


pub struct ReaderOfBams {
    readers: Vec<bam::IndexedReader>,
    filepaths: Vec<String>,
    cur_chr: i32,    
    pub chr_lens: Vec<usize>,
    pub chr_names: Vec<String>,
    prev_reads: Vec<Option<Record>>, // = Vec::with_capacity( n_bams);
    pub record: Option<Record>,
    pub prev_tids: Vec<i32>,
    pub prev_poss: Vec<usize>,
}


impl ReaderOfBams {

    pub fn new(capacity : usize) -> ReaderOfBams { 

        ReaderOfBams {
            readers: Vec::with_capacity( capacity),
            filepaths: Vec::with_capacity( capacity),
            cur_chr : -1,
            chr_lens : Vec::with_capacity( capacity),
            chr_names : Vec::with_capacity( capacity),
            prev_reads : vec![ None; capacity], 
            record : Some(Record::new()),  
            prev_tids : vec![ -1i32; capacity],
            prev_poss : vec![ 0usize; capacity],
        }
    }

    pub fn add_file(&mut self, path: &str) {

    	let bam_metadata = fs::metadata(&path).unwrap_or_else(
    		|_| error!("Cannot read BAM file {}.", path));

    	let index_path = format!("{}.bai", path);
    	let index_metadata = fs::metadata(&index_path).unwrap_or_else(
    		|_| error!("Cannot read BAM index file {}.", index_path));

    	let bam_modified = bam_metadata.modified().unwrap();
    	let index_modified = index_metadata.modified().unwrap();
    	if bam_modified > index_modified {
    		error!("Index file {} is older than its corresponding BAM file.", index_path);
    	}

        self.readers.push(bam::IndexedReader::from_path(path).unwrap_or_else(
        	|_| error!("Cannot open indexed BAM file {}.", path)));
        self.check_header();
        self.filepaths.push(path.to_string());
    }


    fn check_header(&mut self) -> bool {
    	let header_view = self.readers.last()
    		.unwrap_or_else(|| error!("No BAM files have been added"))
    		.header().clone();

        for tid in 0..header_view.target_names().len() {

            let name = String::from_utf8( header_view.target_names()[ tid].to_vec())
            	.unwrap_or_else(|_| error!("Failed to read target names."));
            let chr_len = header_view.target_len( tid as u32)
            	.unwrap_or_else(|| error!("Failed to read target length.")) as usize;

            if self.readers.len() == 1 {
                // Collect names and lengths from first bam file
                self.chr_names.push(name);
                self.chr_lens.push(chr_len);
            } else {
                // Check names and lengths against the first BAM file added
                if name != self.chr_names[tid] || chr_len != self.chr_lens[tid] {
                    // Sanity check
                    error!("Input BAM files have mismatched reference sequences. ('{}')", name);
                }                
            }
        }
        return true;
    }

    // Return chr_index for a given chr_name
    pub fn index_for_name( &self, name : &str) -> Result<usize, bool> {

        let n_regions = self.chr_names.len();

        for i in 0..n_regions {
            //eprintln!( "IND: {},{},{}", &self.chr_names[ i], name, &self.chr_names[ i] == name);
            if &self.chr_names[ i] == name { return Ok( i); }
        }
        return Err(false);
    }

    // Load reads for a region when using bam index files
    // NOTE: It would be possible to limit loading to a window 
    // inside region.
    pub fn load_chr( &mut self, chr : usize) {
        eprintln!("Loading records in region '{}' ({})", self.chr_names[ chr], chr);
        for b in 0..self.readers.len() {
            // Fetch whole chromosome (region)
            self.readers[b].fetch(chr as u32, 0, self.chr_lens[chr] as u64);
        }  

        self.cur_chr = chr as i32;
    }

    // Loads a new record from the bam file from a region and window into self.record
    // Returns true if record was successfully loaded
    // Returns false if there are no more records in the given chr and window
    // In iterative mode (no bam index files) it is not possible to go back from 
    // larger chr indices to smaller ones
    pub fn read_next_record( &mut self, file_index: usize, chr: i32, min_pos : usize, max_pos : usize ) -> bool {

        if chr != self.cur_chr { self.load_chr( chr as usize); }

        let mut tid = -1i32;
        let mut pos = 0usize;
        if chr == -1 {
        	error!("Reader does not support reading unaligned records (chr:-1)." );
        }

        // Use previously saved record?
        if self.prev_reads[ file_index].is_some() {

            tid = self.prev_reads[ file_index].as_ref().unwrap().tid();
            pos = self.prev_reads[ file_index].as_ref().unwrap().pos() as usize;

            if tid == chr && pos <= max_pos && pos >= min_pos {
                //Takes the value out of the option, leaving a None in its place.
                self.record = self.prev_reads[ file_index].take();
                //self.next_pos[ file_index] = 0;
                //self.next_chr[ file_index] = -1;                
                return true;  
            }       
            else if tid > chr || (tid == chr && pos > max_pos) {
                // Saved read is from another (upcoming) region, 
                // Keep it saved
                return false;
            }
            else {
                // Requested chr does not match saved record, discard saved record and
                // continue reading
                eprintln!("WARNING: Sequential read ignored. bam {}", file_index ); 
                //eprintln!( "         Prev: {} pos:{}", tid, pos);
                //eprintln!( "         Cur : {} pos:{}", chr, min_pos);
                self.record = self.prev_reads[ file_index].take();
                tid = -1i32; pos = 0usize;         
            }                         
        }    

        // Read a new record from the bam file
        // until a read from the correct region is found
        // or the end of the file is reached
        while tid < chr || min_pos > pos {

            //eprintln!( "READINGG!");
            match self.readers[file_index].read(&mut self.record.as_mut().unwrap()) {
                Err(bam::Error::TruncatedRecord) => {
                    let bam_path = &self.filepaths[ file_index];
                    error!("Input file '{}' is truncated.", bam_path);
                },
                Err(bam::Error::InvalidRecord) => {
                    let bam_path = &self.filepaths[ file_index];
                    error!("Input file '{}' is corrupted.", bam_path);
                },
                Ok(false) => {                     
                    return false; 
                },
                _ => {}
            }

            tid = self.record.as_ref().unwrap().tid();            
            if tid == -1 {
                // End of the file can contain a lot of unaligned reads
                // that have tid == -1. Check if reading the last region
                // and last window to exit in iterative mode
                if (chr as usize) == self.chr_names.len()-1 && max_pos >= self.chr_lens[ self.chr_lens.len()-1] { return false; }
                continue; // Skip unaligned reads
            } 
            pos = self.record.as_ref().unwrap().pos() as usize;

            // Check positional sorting
            let prev_tid = self.prev_tids[ file_index];
            let prev_pos = self.prev_poss[ file_index];
            
            if prev_tid != -1 && (tid < prev_tid || (tid == prev_tid && prev_pos > pos)) {                 
                error!("BAM file '{}' records are not sorted by position. (prev_tid:{},tid:{},prev_pos:{},pos:{})", self.filepaths[ file_index],prev_tid,tid,prev_pos,pos); 
            }

            self.prev_tids[ file_index] = tid;
            self.prev_poss[ file_index] = pos;

            // Save for later and stop reading
            if tid > chr || (tid == chr && pos > max_pos) { 
                //eprintln!( "SAVING record bam {}: tid{} -> {} ", file_index, tid, pos );
                self.prev_reads[ file_index] = self.record.take(); 
                self.record = Some( Record::new());  // Allocate new record
                return false;
            }                           
        }

        return true;
    }    

    // Returns (-1,0) if no reads have been saved
    pub fn next_read_pos( &self ) -> (i32, usize) {

        let len = self.readers.len();

        let mut min_chr = ::std::i32::MAX;
        let mut min_pos = ::std::usize::MAX;

        for i in 0..len {

            if !self.prev_reads[ i].is_some() { continue; }
            let ref prev_read = self.prev_reads[ i].as_ref().unwrap();

            let tid = prev_read.tid(); 
            if tid == -1 { continue; }

            min_chr = ::std::cmp::min( min_chr, tid );

            if tid <= min_chr { 
                let pos = prev_read.pos() as usize; 
                if tid == min_chr { min_pos = ::std::cmp::min( min_pos, pos ); }
                else { min_pos = pos; }
            }            
        }

        if min_chr == ::std::i32::MAX { return (-1,0); }
        return (min_chr, min_pos);
    }

}
