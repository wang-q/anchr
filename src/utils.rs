use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;

pub fn read_fasta(input: &str) -> anyhow::Result<BTreeMap<String, String>> {
    let mut reader = SeqReader::new(input)?;
    let mut rec = SeqRecord::new();
    let mut seq_of = BTreeMap::new();
    while reader.read_record(&mut rec)? {
        if rec.sequence().is_empty() {
            continue;
        }

        let name = rec.name().to_string();
        let seq = String::from_utf8(rec.sequence().to_vec().to_ascii_uppercase()).unwrap();

        seq_of.insert(name, seq);
    }

    Ok(seq_of)
}

/// Write lines to `outname`, one per line (gzip supported via `pgr::libs::io::writer`).
pub fn write_lines(outname: &str, lines: &[&str]) -> anyhow::Result<()> {
    let mut writer = pgr::libs::io::writer(outname)?;
    for line in lines {
        writeln!(writer, "{}", line)?;
    }
    Ok(())
}

pub fn ucfirst(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

pub fn file_exists(dir: &Path, filename: &str) -> bool {
    let path = dir.join(filename);

    // Check if the file exists
    path.exists()
}

pub fn pause() {
    let mut stdin = std::io::stdin();
    let mut stdout = std::io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}

#[cfg(test)]
mod tests_str {
    use super::*;

    #[test]
    fn test_ucfirst_normal() {
        assert_eq!(ucfirst("hello"), "Hello");
        assert_eq!(ucfirst("rust programming"), "Rust programming");
        assert_eq!(ucfirst("Hello World"), "Hello World");
    }

    #[test]
    fn test_ucfirst_with_numbers() {
        assert_eq!(ucfirst("123abc"), "123abc");
    }

    #[test]
    fn test_ucfirst_empty() {
        assert_eq!(ucfirst(""), "");
    }

    #[test]
    fn test_ucfirst_special_characters() {
        assert_eq!(ucfirst("!hello"), "!hello");
        assert_eq!(ucfirst("@rust"), "@rust");
    }

    #[test]
    fn test_ucfirst_unicode() {
        assert_eq!(ucfirst("你好"), "你好");
        assert_eq!(ucfirst("こんにちは"), "こんにちは"); // 日文字符
    }

    #[test]
    fn test_ucfirst_leading_spaces() {
        assert_eq!(ucfirst("   leading spaces"), "   leading spaces");
    }

    #[test]
    fn test_ucfirst_single_character() {
        assert_eq!(ucfirst("a"), "A");
        assert_eq!(ucfirst("Z"), "Z");
    }
}

#[cfg(test)]
mod tests_file {
    #[test]
    fn test_file_exists() {
        // Create a temporary directory
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path().to_path_buf();

        // Test case for an existing file
        let filename = "test_file.txt";
        let file_path = dir_path.join(filename);

        // Create a test file
        std::fs::write(&file_path, "Hello, World!").expect("Unable to write to file");

        // Check if the file exists
        assert!(crate::file_exists(&dir_path, filename));
    }

    #[test]
    fn test_file_does_not_exist() {
        // Create a temporary directory
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path().to_path_buf();

        // Check for a non-existent file
        let filename = "non_existent_file.txt";
        assert!(!crate::file_exists(&dir_path, filename));
    }

    #[test]
    fn test_empty_directory() {
        // Create a temporary directory
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path().to_path_buf();

        // Check for a non-existent file in an empty directory
        let filename = "empty_file.txt";
        assert!(!crate::file_exists(&dir_path, filename));
    }

    #[test]
    fn test_file_with_special_characters() {
        // Create a temporary directory
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path().to_path_buf();

        // Create a file with special characters in its name
        let filename = "special_file_@#$.txt";
        let file_path = dir_path.join(filename);
        std::fs::write(&file_path, "Special characters!").expect("Unable to write to file");

        // Check if the file exists
        assert!(crate::file_exists(&dir_path, filename));
    }
}
