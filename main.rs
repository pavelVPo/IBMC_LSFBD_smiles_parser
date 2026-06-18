mod ps_supps;                                        //
//use std::fmt;                                     //   Output
use wasm_bindgen::prelude::*;                      //    JavaScript bindings
use ps_supps::*;                                  //
//_______________________________________________//

// #red *, anything -> any atom only will be considered

fn main() {

    // Test
    let mut this_structure = Structure::default();
    let smiles: String = "]CCCC".to_string();

    this_structure = parse_smiles(&smiles, this_structure);
    
    println!("{:?}", this_structure);
}