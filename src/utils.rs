use lazy_static::lazy_static;
use std::collections::HashMap;
use std::io::stdin;

lazy_static! {
    static ref CODON_TABLE: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        m.insert("F", "TTT");
        m.insert("L", "CTG");
        m.insert("S", "AGC");
        m.insert("Y", "TAC");
        m.insert("*", "TGA");
        m.insert("C", "TGC");
        m.insert("W", "TGG");
        m.insert("P", "CCG");
        m.insert("H", "CAC");
        m.insert("Q", "CAG");
        m.insert("R", "CGG");
        m.insert("I", "ATC");
        m.insert("M", "ATG");
        m.insert("T", "ACC");
        m.insert("N", "AAC");
        m.insert("K", "AAG");
        m.insert("V", "GTG");
        m.insert("A", "GCC");
        m.insert("D", "GAC");
        m.insert("E", "GAG");
        m.insert("G", "GGC");
        m
    };
}

pub fn print_translation(translation: &str) {
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

pub fn read_mutations(limit: u32, seq_len: usize) -> Vec<(u32, String)> {
    let mut changes = Vec::new();
    println!(
        "Enter the position, space, and the amino acid to mutate to (or type 'end' to finish)"
    );

    let mut end_flag = false; // Flag to indicate if "end" was typed

    for i in 0..limit {
        if end_flag {
            break; // Break out of the outer loop if "end" was typed
        }

        loop {
            let mut input = String::new();
            println!("Mutation {}:", i + 1);
            stdin().read_line(&mut input).expect("Failed to read line");

            let input = input.trim();
            if input == "end" {
                end_flag = true; // Set the flag to true
                break; // Break out of the inner loop
            }
            let input: Vec<&str> = input.split_whitespace().collect();
            if input.len() != 2 {
                println!("Invalid input format. Please enter in the format 'num AA'.");
                continue;
            }
            if let Ok(num) = input[0].parse::<u32>() {
                if num == 0 || num > seq_len as u32 {
                    println!(
                        "Position out of range. Please enter a position between 1 and {}.",
                        seq_len
                    );
                    continue;
                }
                let aa = input[1].to_string();
                if !CODON_TABLE.contains_key(aa.as_str()) {
                    println!("Invalid amino acid. Please enter a valid amino acid.");
                    continue;
                }
                changes.push((num - 1, aa));
                break; // Break out of the inner loop to proceed to the next mutation
            } else {
                println!("Invalid number format. Please enter a valid number.");
            }
        }
    }

    changes
}

pub fn design_primers(
    changes: Vec<(u32, String)>,
    sequence: &str,
) -> Result<(String, String), String> {
    let mut new_sequence = sequence.to_string();

    for (num, aa) in &changes {
        if let Some(codon) = CODON_TABLE.get(aa.as_str()) {
            let start = (num * 3) as usize;
            let end = start + 3;
            if end <= sequence.len() {
                new_sequence.replace_range(start..end, codon);
            }
        }
    }

    let mut start = (changes.iter().map(|(num, _)| num).min().unwrap() * 3) as usize;
    let mut end = (changes.iter().map(|(num, _)| num).max().unwrap() * 3 + 3) as usize;
    println!("{}", start);
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
        } else if end < new_sequence.len() && !start_turn {
            end += 1;
            start_turn = true;
        } else {
            return Err("Desired mutation(s) is too close to the start or end. Please input new sequence with the intended mutation(s) centered.".to_string());
        }
    }

    let primer = new_sequence[start..end].to_string();
    let rev_comp = reverse_transcribe(&primer);

    Ok((rev_comp, primer))
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
