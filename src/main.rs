mod utils;

use std::io::stdin;

fn main() {
    let mut sequence = String::new();
    println!("Input DNA Sequence 35 bp or larger");
    stdin()
        .read_line(&mut sequence)
        .expect("Failed to read line");

    let sequence = sequence.trim().to_string();

    // Check if the sequence is less than 35 bp
    if sequence.len() < 35 {
        eprintln!("Error: The input DNA sequence must be at least 35 base pairs long.");
        return;
    }

    let translation = utils::translate(&sequence);
    utils::print_translation(&translation);

    let limit = 5;
    let mutations = utils::read_mutations(limit, translation.len());

    match utils::design_primers(mutations, &sequence) {
        Ok((top, btm)) => println!("Top: {}\nBottom: {}", top, btm),
        Err(e) => println!("Error: {}", e),
    }
}
