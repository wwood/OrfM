extern crate assert_cli;

#[cfg(test)]
mod tests {
    #[test]
    fn test_inseqs_fna() {
        let expected = std::fs::read_to_string("tests/data/inseqs.fna.expected").unwrap();
        assert_cli::Assert::main_binary()
            .with_args(&["-m", "72", "-c", "4", "tests/data/inseqs.fna"])
            .stdout()
            .is(expected.trim_end())
            .unwrap();
    }

    #[test]
    fn test_inseqs_fna_process_substitution() {
        let expected = std::fs::read_to_string("tests/data/inseqs.fna.expected").unwrap();
        let bin = env!("CARGO_BIN_EXE_orfm");
        let output = std::process::Command::new("bash")
            .arg("-c")
            .arg(format!("{} -m 72 -c 4 <(cat tests/data/inseqs.fna)", bin))
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "orfm failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert_eq!(stdout.trim_end(), expected.trim_end());
    }
}
