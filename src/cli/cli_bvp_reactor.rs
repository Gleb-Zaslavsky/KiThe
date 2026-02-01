use crate::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use crate::ReactorsBVP::task_parser_reactor_BVP::create_template;
use crate::cli::reactor_help::{REACTOR_ENG_HELPER, REACTOR_RU_HELPER};
use RustedSciThe::Utils::task_parser::pretty_print_map;
use dialoguer::{Confirm, Input, Select, theme::ColorfulTheme};
use std::path::PathBuf;

pub fn reactor_menu() {
    let theme = ColorfulTheme::default();
    loop {
        let items = vec![
            "Solve from file",
            "Auto-discover problem files",
            "Generate template",
            "Read help (eng)",
            "Read help (ru)",
            "Back to main menu",
        ];

        let selection = Select::with_theme(&theme)
            .with_prompt("=== Reactor BVP Problems ===")
            .items(&items)
            .default(0)
            .interact()
            .unwrap();

        match selection {
            0 => solve_from_file(),
            1 => auto_solve_problems(),
            2 => {
                create_template();
                println!("Template generated successfully!");
            }
            3 => show_help_english(),
            4 => show_help_russian(),
            5 => break,
            _ => {}
        }
    }
}

fn solve_from_file() {
    let theme = ColorfulTheme::default();
    let file_path: String = Input::with_theme(&theme)
        .with_prompt("Enter file path")
        .interact_text()
        .unwrap();
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

pub fn solve_from_file_dialog(path: std::path::PathBuf) {
    let theme = ColorfulTheme::default();
    let mut reactor = SimpleReactorTask::new();
    // Test parsing
    let mut parser = reactor
        .parse_file(Some(path))
        .expect("Failed to parse file");
    match parser.parse_document() {
        Ok(result) => {
            println!("Document parsed successfully");
            pretty_print_map(&result);

            let start = Confirm::with_theme(&theme)
                .with_prompt("Start calculation?")
                .default(true)
                .interact()
                .unwrap();

            if start {
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
    println!("{}", REACTOR_ENG_HELPER);
    let _ = std::io::stdin().read_line(&mut String::new());
}

fn show_help_russian() {
    println!("\n=== Справка по реакторам (Русский) ===");
    println!("{}", REACTOR_RU_HELPER);
    let _ = std::io::stdin().read_line(&mut String::new());
}
