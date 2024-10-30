use std::collections::HashMap;
use std::io::stdin;

pub fn print_translation(translation: String) {
    let mut printed: String = String::new();
    for (i, c) in translation.chars().enumerate() {
        if i % 10 == 0 && i > 0 {
            printed.push_str("\n");
        }
        printed.push_str(&format!("{}: {}\t", i + 1, c));
    }
    println!("{}", printed);
}

pub fn translate(sequence: &str) -> String {
    let codon_table: HashMap<&str, &str> = [
        ("TTT", "F"),
        ("TTC", "F"),
        ("TTA", "L"),
        ("TTG", "L"),
        ("TCT", "S"),
        ("TCC", "S"),
        ("TCA", "S"),
        ("TCG", "S"),
        ("TAT", "Y"),
        ("TAC", "Y"),
        ("TAA", "*"),
        ("TAG", "*"),
        ("TGT", "C"),
        ("TGC", "C"),
        ("TGA", "*"),
        ("TGG", "W"),
        ("CTT", "L"),
        ("CTC", "L"),
        ("CTA", "L"),
        ("CTG", "L"),
        ("CCT", "P"),
        ("CCC", "P"),
        ("CCA", "P"),
        ("CCG", "P"),
        ("CAT", "H"),
        ("CAC", "H"),
        ("CAA", "Q"),
        ("CAG", "Q"),
        ("CGT", "R"),
        ("CGC", "R"),
        ("CGA", "R"),
        ("CGG", "R"),
        ("ATT", "I"),
        ("ATC", "I"),
        ("ATA", "I"),
        ("ATG", "M"),
        ("ACT", "T"),
        ("ACC", "T"),
        ("ACA", "T"),
        ("ACG", "T"),
        ("AAT", "N"),
        ("AAC", "N"),
        ("AAA", "K"),
        ("AAG", "K"),
        ("AGT", "S"),
        ("AGC", "S"),
        ("AGA", "R"),
        ("AGG", "R"),
        ("GTT", "V"),
        ("GTC", "V"),
        ("GTA", "V"),
        ("GTG", "V"),
        ("GCT", "A"),
        ("GCC", "A"),
        ("GCA", "A"),
        ("GCG", "A"),
        ("GAT", "D"),
        ("GAC", "D"),
        ("GAA", "E"),
        ("GAG", "E"),
        ("GGT", "G"),
        ("GGC", "G"),
        ("GGA", "G"),
        ("GGG", "G"),
    ]
    .iter()
    .cloned()
    .collect();

    let mut translation: String = String::new();

    for i in (0..sequence.len()).step_by(3) {
        if let Some(codon) = sequence.get(i..i + 3) {
            if let Some(amino_acid) = codon_table.get(codon) {
                translation.push_str(amino_acid);
            }
        }
    }

    translation
}

pub fn read_mutations() -> Vec<(u32, String)> {
    let mut changes = Vec::new();
    println!(
        "Enter the position, space, and the amino acid to mutate to (or type 'end' to finish)"
    );

    let mut i = 1;
    loop {
        let mut input = String::new();
        println!("Mutation {}:", i);
        stdin().read_line(&mut input).expect("Failed to read line");

        let input = input.trim();
        if input == "end" {
            break;
        }
        let input: Vec<&str> = input.split_whitespace().collect();
        if let Ok(num) = input[0].parse::<u32>() {
            let aa = input[1].to_string();
            changes.push((num - 1, aa));
            i += 1;
        }
    }

    changes
}

pub fn design_primers(changes: Vec<(u32, String)>, sequence: &str) -> (String, String) {
    let codon_table: HashMap<&str, &str> = [
        ("F", "TTT"),
        ("L", "CTG"),
        ("S", "AGC"),
        ("Y", "TAC"),
        ("*", "TGA"),
        ("C", "TGC"),
        ("W", "TGG"),
        ("P", "CCG"),
        ("H", "CAC"),
        ("Q", "CAG"),
        ("R", "CGG"),
        ("I", "ATC"),
        ("M", "ATG"),
        ("T", "ACC"),
        ("N", "AAC"),
        ("K", "AAG"),
        ("V", "GTG"),
        ("A", "GCC"),
        ("D", "GAC"),
        ("E", "GAG"),
        ("G", "GGC"),
    ]
    .iter()
    .cloned()
    .collect();

    let mut new_sequence = sequence.to_string();

    for (num, aa) in &changes {
        if let Some(codon) = codon_table.get(aa.as_str()) {
            let start = (num * 3) as usize;
            let end = start + 3;
            if end <= sequence.len() {
                new_sequence.replace_range(start..end, codon);
            }
        }
    }

    let mut start = (changes.iter().map(|(num, _)| num).min().unwrap() * 3) as usize;
    let mut end = (changes.iter().map(|(num, _)| num).max().unwrap() * 3 + 3) as usize;

    let mut start_turn = true;

    while end - start < 25 || {
        let pmis = percent_mismatch(&sequence[start..end], &new_sequence[start..end]);
        let pgc = percent_gc(&new_sequence[start..end]);
        let tm = 81.5 + 0.41 * pgc - 675.0 / (end - start) as f32 - pmis;
        tm < 78.0
    } {
        if start > 0 && start_turn {
            start -= 1;
            start_turn = false;
        } else if end < new_sequence.len() {
            end += 1;
            start_turn = true;
        }
    }

    let primer = new_sequence[start..end].to_string();
    let rev_comp = reverse_transcribe(&primer);

    (rev_comp, primer)
}

fn percent_mismatch(old_seq: &str, new_seq: &str) -> f32 {
    let mut mismatch = 0;
    for (old_char, new_char) in old_seq.chars().zip(new_seq.chars()) {
        if old_char != new_char {
            mismatch += 1;
        }
    }
    (mismatch as f32 / old_seq.len() as f32) * 100.0
}

fn percent_gc(sequence: &str) -> f32 {
    let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
    let total_count = sequence.len();
    if total_count == 0 {
        return 0.0;
    }
    (gc_count as f32 / total_count as f32) * 100.0
}

fn reverse_transcribe(sequence: &str) -> String {
    let key: HashMap<char, char> = [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]
        .iter()
        .cloned()
        .collect();

    let mut rev_comp = String::new();

    for nuc in sequence.chars().rev() {
        if let Some(&comp) = key.get(&nuc) {
            rev_comp.push(comp);
        }
    }
    rev_comp
}
