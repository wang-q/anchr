use bio::io::fasta;
use std::collections::BTreeMap;

pub fn read_fasta(input: &str) -> BTreeMap<String, String> {
    let mut reader = intspan::reader(input);
    let fa_in = fasta::Reader::new(reader);

    let mut seq_of = BTreeMap::new();
    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result.unwrap();

        if record.is_empty() {
            continue;
        }

        let name = record.id().to_string();
        let seq = String::from_utf8(record.seq().to_vec().to_ascii_uppercase()).unwrap();

        seq_of.insert(name, seq);
    }

    seq_of
}
