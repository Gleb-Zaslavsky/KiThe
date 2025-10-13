fn main() {
    #[cfg(windows)]
    {
        if std::path::Path::new("src/assets/icon.ico").exists() {
            let mut res = winres::WindowsResource::new();
            res.set_icon("src/assets/icon.ico");
            if let Err(e) = res.compile() {
                println!("cargo:warning=Failed to compile icon resource: {}", e);
            }
        } else {
            println!(
                "cargo:warning=Icon file 'src/assets/icon.ico' not found. Skipping icon embedding."
            );
        }
    }
}
