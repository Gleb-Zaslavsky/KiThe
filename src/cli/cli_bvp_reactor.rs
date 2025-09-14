use crate::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use crate::ReactorsBVP::task_parser_reactor_BVP::create_template;
use crate::cli::reactor_help::{REACTOR_ENG_HELPER, REACTOR_RU_HELPER};
use RustedSciThe::Utils::task_parser::pretty_print_map;
use std::io::{self, Write};
use std::path::PathBuf;

pub fn reactor_menu() {
    loop {
        println!("\n=== Reactor BVP Problems ===");
        println!("\x1b[33m1. Solve from file\x1b[0m");
        println!("\x1b[33m2. Auto-discover problem files\x1b[0m");
        println!("\x1b[33m3. Generate template\x1b[0m");
        println!("\x1b[33m4. Read help (eng)\x1b[0m");
        println!("\x1b[33m5. Read help (ru)\x1b[0m");
        println!("\x1b[33m0. Back to main menu\x1b[0m");
        print!("\x1b[36mEnter your choice: \x1b[0m");
        io::stdout().flush().unwrap();

        let choice = get_user_input();
        match choice.trim() {
            "1" => solve_from_file(),
            "2" => auto_solve_problems(),
            "3" => {
                create_template();
                println!("Template generated successfully!");
            }
            "4" => show_help_english(),
            "5" => show_help_russian(),
            "0" => break,
            _ => println!("Invalid choice. Please try again."),
        }
    }
}

fn solve_from_file() {
    print!("\x1b[36mEnter file path: \x1b[0m");
    io::stdout().flush().unwrap();
    let file_path = get_user_input();
    let path = PathBuf::from(file_path.trim());

    if path.exists() {
        solve_from_file_dialog(path);
    } else {
        println!("File not found: {}", file_path.trim());
    }
}

fn auto_solve_problems() {
    use std::{env, fs};

    let current_dir = env::current_dir().expect("Failed to get current directory");
    println!("Searching for problem files in: {:?}", current_dir);

    let mut found_files = false;

    if let Ok(entries) = fs::read_dir(&current_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                let path = entry.path();
                if let Some(filename) = path.file_name() {
                    let filename_str = filename.to_string_lossy();
                    if filename_str.starts_with("problem") && path.is_file() {
                        println!("Found problem file: {:?}", path);
                        solve_from_file_dialog(path);
                        found_files = true;
                    }
                }
            }
        }
    }

    if !found_files {
        println!("No files starting with 'problem' found in current directory.");
    }
}

fn get_user_input() -> String {
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read input");
    input
}

pub fn solve_from_file_dialog(path: std::path::PathBuf) {
    let mut reactor = SimpleReactorTask::new();
    // Test parsing
    let mut parser = reactor
        .parse_file(Some(path))
        .expect("Failed to parse file");
    match parser.parse_document() {
        Ok(result) => {
            println!("Document parsed successfully");
            pretty_print_map(&result);

            print!("\x1b[36mStart calculation? (y/n): \x1b[0m");
            io::stdout().flush().unwrap();
            let choice = get_user_input();

            if choice.trim().to_lowercase() == "y" || choice.trim().to_lowercase() == "yes" {
                reactor.solve_from_map(parser);
            } else {
                println!("Calculation cancelled. Returning to menu.");
            }
        }
        Err(_) => {
            let err = parser.get_error().unwrap();
            println!("Error parsing document: {}", err);
        }
    }
}

fn show_help_english() {
    println!("\n=== Reactor BVP Help (English) ===");
    println!("\nPress Enter to return to menu...");
    println!("{}", REACTOR_ENG_HELPER);
    let _ = get_user_input();
}

fn show_help_russian() {
    println!("\n=== Справка по реакторам (Русский) ===");
    println!("\nНажмите Enter для возврата в меню...");
    println!("{}", REACTOR_RU_HELPER);
    let _ = get_user_input();
}
