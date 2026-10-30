use serde::Deserialize;
use std::{env, fs, path::PathBuf};
use tabled::{Table, Tabled};

#[derive(Deserialize)]
struct Campaign {
    run_id: String,
    suite: String,
    created_utc: String,
    duration_seconds: f64,
    git_revision: String,
    working_tree_dirty: bool,
    host: String,
    rustc: String,
    cargo: String,
    build_jobs: usize,
    tcc_available: bool,
    steps: Vec<Step>,
    criterion_rows: Vec<CriterionRow>,
    diagnostic_tables_markdown: String,
}

#[derive(Deserialize, Tabled)]
struct Step {
    #[tabled(rename = "Step")]
    id: String,
    #[tabled(rename = "Suite")]
    suite: String,
    #[tabled(rename = "Status")]
    status: String,
    #[tabled(rename = "Seconds")]
    duration_seconds: f64,
    #[tabled(rename = "Exit")]
    #[tabled(display("tabled::derive::display::option", "-"))]
    exit_code: Option<i32>,
    #[tabled(rename = "Technical log")]
    technical_log: String,
    #[serde(skip)]
    #[tabled(skip)]
    description: String,
    #[serde(skip)]
    #[tabled(skip)]
    reason: String,
    #[serde(skip)]
    #[tabled(skip)]
    command: String,
}

#[derive(Deserialize, Tabled)]
struct CriterionRow {
    #[tabled(rename = "Benchmark")]
    #[serde(default)]
    benchmark: String,
    #[tabled(rename = "Scenario")]
    fixture: String,
    #[tabled(rename = "Solver/frontend/layout")]
    route: String,
    #[tabled(rename = "Median ms")]
    median_ms: f64,
    #[tabled(rename = "95% low ms")]
    ci_low_ms: f64,
    #[tabled(rename = "95% high ms")]
    ci_high_ms: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args_os().skip(1);
    let input = PathBuf::from(args.next().ok_or("missing campaign.json path")?);
    let output = PathBuf::from(args.next().ok_or("missing summary.md path")?);
    let bytes = fs::read(input)?;
    let json = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(&bytes);
    let campaign: Campaign = serde_json::from_slice(json)?;
    let failed = campaign
        .steps
        .iter()
        .filter(|step| step.status == "failed" || step.status == "interrupted")
        .count();

    let mut report = format!(
        "# ReactorBVP release campaign\n\n\
         - Run: `{}` (`{}`)\n\
         - Created (UTC): `{}`\n\
         - Duration: `{:.3} s`\n\
         - Revision: `{}` (working tree dirty: `{}`)\n\
         - Host/toolchain: `{}` / `{}` / `{}`\n\
         - Build jobs: `{}`; TCC available: `{}`\n\
         - Steps: {}; failed/interrupted: `{}`\n\n\
         ## Campaign status\n\n```text\n{}\n```\n",
        campaign.run_id,
        campaign.suite,
        campaign.created_utc,
        campaign.duration_seconds,
        campaign.git_revision,
        campaign.working_tree_dirty,
        campaign.host,
        campaign.rustc,
        campaign.cargo,
        campaign.build_jobs,
        campaign.tcc_available,
        campaign.steps.len(),
        failed,
        Table::new(&campaign.steps)
    );

    report.push_str(
        "\n## Criterion release wall time\n\nMedian and 95% estimate interval, milliseconds.\n\n",
    );
    if campaign.criterion_rows.is_empty() {
        report
            .push_str("No Criterion estimates were found; inspect the technical benchmark log.\n");
    } else {
        report.push_str(&Table::new(&campaign.criterion_rows).to_string());
        report.push('\n');
    }

    report.push_str("\n## Detailed stage and native telemetry tables\n\n");
    if campaign.diagnostic_tables_markdown.trim().is_empty() {
        report.push_str("No Tabled output was extracted; inspect the raw benchmark log.\n");
    } else {
        report.push_str(&campaign.diagnostic_tables_markdown);
    }

    let details = campaign
        .steps
        .iter()
        .map(|step| {
            let reason = if step.reason.is_empty() {
                String::new()
            } else {
                format!("\n  Result: {}", step.reason)
            };
            format!(
                "- **{}:** {}{}\n  Command: `{}`",
                step.id, step.description, reason, step.command
            )
        })
        .collect::<Vec<_>>();
    if !details.is_empty() {
        report.push_str("\n## Step details\n\n");
        report.push_str(&details.join("\n"));
        report.push('\n');
    }
    fs::write(output, report)?;
    Ok(())
}
