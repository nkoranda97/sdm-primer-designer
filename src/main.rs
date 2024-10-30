mod utils;

use std::io::stdin;

fn main() {
    let mut sequence = String::new();
    println!("Input DNA Sequence");
    stdin()
        .read_line(&mut sequence)
        .expect("Failed to read line");

    let translation = utils::translate(&sequence);
    utils::print_translation(translation);

    let mutations = utils::read_mutations();

    let (btm, top) = utils::design_primers(mutations, &sequence);

    println!("top: {}\nbtm: {}", top, btm);
}
