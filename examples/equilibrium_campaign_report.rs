use serde::Deserialize;
use std::{env, fs, path::PathBuf};
use tabled::{Table, Tabled};

#[derive(Deserialize)]
struct Campaign {
    run_id: String,
    suite: String,
    profile: String,
    created_utc: String,
    duration_seconds: f64,
    git_revision: String,
    working_tree_dirty: bool,
    host: String,
    rustc: String,
    cargo: String,
    build_jobs: usize,
    steps: Vec<Step>,
}

#[derive(Deserialize, Tabled)]
struct Step {
    #[tabled(rename = "Step")]
    id: String,
    #[tabled(rename = "Profile")]
    profile: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Seconds")]
    duration_seconds: f64,
    #[tabled(rename = "Exit")]
    #[tabled(display("tabled::derive::display::option", "-"))]
    exit_code: Option<i32>,
    #[tabled(rename = "Workload")]
    #[tabled(display = "compact_workload")]
    workload: String,
    #[tabled(rename = "Log")]
    technical_log: String,
    #[serde(skip)]
    #[tabled(skip)]
    reason: String,
}

fn compact_workload(workload: &String) -> String {
    const MAX_CHARS: usize = 34;
    if workload.chars().count() <= MAX_CHARS {
        workload.clone()
    } else {
        let prefix = workload.chars().take(MAX_CHARS - 3).collect::<String>();
        format!("{prefix}...")
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args_os().skip(1);
    let input = PathBuf::from(args.next().ok_or("missing campaign.json path")?);
    let output = PathBuf::from(args.next().ok_or("missing summary.md path")?);
    let bytes = fs::read(input)?;
    let json = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(&bytes);
    let campaign: Campaign = serde_json::from_slice(json)?;
    let table = Table::new(&campaign.steps).to_string();
    let failed = campaign
        .steps
        .iter()
        .filter(|step| step.status != "passed")
        .count();
    let details = campaign
        .steps
        .iter()
        .filter(|step| !step.reason.is_empty())
        .map(|step| format!("- **{}:** {}", step.id, step.reason))
        .collect::<Vec<_>>();
    let report = format!(
        "# Equilibrium campaign\n\n\
         - Run: `{}`\n\
         - Suite/profile: `{}` / `{}`\n\
         - Created (UTC): `{}`\n\
         - Duration: `{:.3} s`\n\
         - Revision: `{}` (working tree dirty: `{}`)\n\
         - Host: `{}`\n\
         - Toolchain: `{}`; `{}`\n\
         - Cargo build jobs: `{}`\n\
         - Steps: {}; non-passing: {}\n\n\
         ```text\n{}\n```\n",
        campaign.run_id,
        campaign.suite,
        campaign.profile,
        campaign.created_utc,
        campaign.duration_seconds,
        campaign.git_revision,
        campaign.working_tree_dirty,
        campaign.host,
        campaign.rustc,
        campaign.cargo,
        campaign.build_jobs,
        campaign.steps.len(),
        failed,
        table
    );
    let report = if details.is_empty() {
        report
    } else {
        format!("{report}\n## Step details\n\n{}\n", details.join("\n"))
    };
    fs::write(output, report)?;
    Ok(())
}
