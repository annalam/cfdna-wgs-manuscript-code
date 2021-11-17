use std::str;
use std::fs::File;
use bio::io::fasta;

pub struct RefGenomeReader {
    filepath: String,
    genome_reader: fasta::IndexedReader<File>,
}

impl RefGenomeReader {

    pub fn new( genome_fasta_path: &str) -> RefGenomeReader {
        RefGenomeReader {
            filepath: genome_fasta_path.to_string(),
            genome_reader: fasta::IndexedReader::from_file(&genome_fasta_path).unwrap_or_else(|_| error!("Could not open genome FASTA file '{}'.", &genome_fasta_path))
        }
    }

    pub fn load_chromosome_seq( &mut self, chr_name : &str ) -> Vec<u8> {
        let mut chr_seq: Vec<u8> = Vec::new();
        self.genome_reader.fetch_all(&chr_name);
        self.genome_reader.read(&mut chr_seq);
        eprintln!("Loaded reference chromosome {} (length {}).", &chr_name,
        	chr_seq.len());
        chr_seq
    }

}

