use matrix_mult2::matrix::{Matrix, multiply, strassen, check};
use std::env;
use std::time::Instant;

fn get_time(size: usize, do_print: bool) {
    let a = Matrix::random(size);
    let b = Matrix::random(size);
    let now = Instant::now();
    let res_a = multiply(&a, &b);
    let time_standard = now.elapsed().as_secs_f64();
    let now = Instant::now();
    let res_b = strassen(&a, &b);
    let time_strassen = now.elapsed().as_secs_f64();
    if do_print {
        println!("Matrix a:");
        a.print();
        println!("Matrix b:");
        b.print();
        println!("Result standard:");
        res_a.print();
        println!("Result Strassen:");
        res_b.print();
    }
    assert!(check(&res_a, &res_b));
    println!("Size: {} Time standard: {} Strassen: {}",
        size,
        time_standard,
        time_strassen
    );

}

pub fn do_test() {
    let mut size = 3;
    while size < 20 {
        get_time(size, true);
        size *= 2;
    }
}

fn main() {
    let size;
    let mode;
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        do_test();
        println!("Error: invalid arguments.");
        println!("The first argument selects the program mode. 0 - experiment mode used in the report, 1 - prints matrices, 2 - doesn't print matrices.");
        println!("The second argument selects the order of matrices. It is only checked for modes 1 and 2");
        println!("Example: program_name.exe 1 256");
        return;
    }
    mode = args[1].parse::<usize>().unwrap();

    match mode {
        1 => {
            size = args[2].parse::<usize>().unwrap();
            get_time(size, true);
        },
        2 => {
            size = args[2].parse::<usize>().unwrap();
            get_time(size, false);
        },
        _ => {
            do_test();
            return;
        }
    }
}
