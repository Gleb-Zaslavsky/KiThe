use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn show_image(image_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let image_path = find_image_in_project(image_name)?;
    open_with_default_viewer(&image_path)?;
    Ok(())
}

fn find_image_in_project(image_name: &str) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let current_dir = env::current_dir()?;

    // Find KiThe directory
    let kithe_dir = find_kithe_directory(&current_dir)?;

    // Search for the image file
    find_image_recursive(&kithe_dir, image_name)
}

fn find_kithe_directory(start_path: &Path) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let mut current = start_path;

    loop {
        if current
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.starts_with("KiThe"))
            .unwrap_or(false)
        {
            return Ok(current.to_path_buf());
        }

        match current.parent() {
            Some(parent) => current = parent,
            None => return Err("KiThe directory not found".into()),
        }
    }
}

fn find_image_recursive(
    dir: &Path,
    image_name: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                if filename == image_name
                    || filename.starts_with(&format!(
                        "{}",
                        image_name.trim_end_matches(|c: char| !c.is_alphanumeric())
                    ))
                {
                    return Ok(path);
                }
            }
        } else if path.is_dir() {
            if let Ok(found) = find_image_recursive(&path, image_name) {
                return Ok(found);
            }
        }
    }

    Err(format!("Image '{}' not found", image_name).into())
}

fn open_with_default_viewer(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(target_os = "windows")]
    {
        Command::new("cmd")
            .args(["/C", "start", "", path.to_str().unwrap()])
            .spawn()?;
    }

    #[cfg(target_os = "macos")]
    {
        Command::new("open").arg(path).spawn()?;
    }

    #[cfg(target_os = "linux")]
    {
        Command::new("xdg-open").arg(path).spawn()?;
    }

    Ok(())
}
