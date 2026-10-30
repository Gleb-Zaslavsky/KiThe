use std::error::Error;

use KiThe::Thermodynamics::DBhandlers::library_audit::audit_library_records_with_telemetry;
use KiThe::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use KiThe::Thermodynamics::thermo_lib_api::ThermoData;

fn main() -> Result<(), Box<dyn Error>> {
    let (mut repository, load) = ThermoData::open_with_telemetry()?;
    println!(
        "repository open: cache={:?}, total_ms={:.3}",
        load.cache_status,
        load.total.unwrap_or_default().as_secs_f64() * 1_000.0
    );

    let selected_keys = vec!["CO".to_string(), "CH4".to_string()];
    let (audit, audit_timing) =
        audit_library_records_with_telemetry(&repository, "NASA_gas", &selected_keys)?;
    println!(
        "audit NASA_gas: requested={}, valid={}, failed={}, total_ms={:.3}",
        audit_timing.records_requested,
        audit_timing.records_valid,
        audit_timing.records_failed,
        audit_timing.total.as_secs_f64() * 1_000.0
    );
    audit.write_json("nasa-gas-selected-audit.json")?;

    if audit.is_clean() {
        let mut production =
            SubsData::from_thermo_repository(ThermoData::try_default_repository()?);
        production.set_substances(selected_keys);
        production.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        production.search_substances()?;
        println!("production selection committed transactionally");
    } else {
        println!("production selection skipped: audit found malformed rows");
    }

    repository.close();
    Ok(())
}
