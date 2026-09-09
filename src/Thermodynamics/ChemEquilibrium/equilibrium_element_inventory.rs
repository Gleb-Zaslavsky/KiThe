//! Validated closed elemental inventories for equilibrium problem assembly.
//!
//! A formal carrier such as `C5H6N7: 1` is accepted only as compact syntax for
//! an elemental contribution. It is never a thermochemical species, a phase
//! component, or a numerical trace coordinate. Once constructed, an
//! [`ElementInventory`] is the canonical physical vector `b` in `A^T n = b`.

use std::collections::BTreeMap;

use thiserror::Error;

use crate::Kinetics::molmass::{is_known_element, parse_formula};

/// One formula-shaped contribution to a closed elemental inventory.
///
/// `formula` is input syntax only. Its parsed atoms are accumulated into the
/// inventory and the string is intentionally not retained as a species-like
/// identity by [`ElementInventory`].
#[derive(Debug, Clone, PartialEq)]
pub struct FormalElementCarrier {
    formula: String,
    amount: f64,
}

impl FormalElementCarrier {
    /// Creates one finite, non-negative formal elemental contribution.
    pub fn new(formula: impl Into<String>, amount: f64) -> Result<Self, ElementInventoryError> {
        let formula = formula.into();
        if formula.trim().is_empty() {
            return Err(ElementInventoryError::EmptyFormula);
        }
        validate_amount("formal carrier", amount)?;
        Ok(Self { formula, amount })
    }

    pub fn formula(&self) -> &str {
        &self.formula
    }

    pub fn amount(&self) -> f64 {
        self.amount
    }
}

/// Canonical, validated elemental inventory `b` for a closed equilibrium
/// system.
///
/// Entries are held in a `BTreeMap`, so caller order has no physical meaning.
/// Exact zero entries are omitted after aggregation; every retained amount is
/// finite and strictly positive.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementInventory {
    amounts: BTreeMap<String, f64>,
}

impl ElementInventory {
    /// Builds an inventory directly from element amounts.
    pub fn from_amounts<I, S>(entries: I) -> Result<Self, ElementInventoryError>
    where
        I: IntoIterator<Item = (S, f64)>,
        S: Into<String>,
    {
        let mut amounts = BTreeMap::new();
        for (element, amount) in entries {
            let element = normalize_element(element.into())?;
            validate_amount(&element, amount)?;
            *amounts.entry(element).or_insert(0.0) += amount;
        }
        Self::from_canonical_amounts(amounts)
    }

    /// Parses formal carriers and adds their elemental contributions.
    ///
    /// The established shared formula parser provides syntax handling; this
    /// boundary additionally verifies that every parsed token is a real
    /// periodic-table element before it can enter a physical inventory.
    pub fn from_formal_carriers<I>(carriers: I) -> Result<Self, ElementInventoryError>
    where
        I: IntoIterator<Item = FormalElementCarrier>,
    {
        let mut amounts = BTreeMap::new();
        for carrier in carriers {
            let parsed = parse_formula(carrier.formula.clone(), None).map_err(|error| {
                ElementInventoryError::MalformedFormula {
                    formula: carrier.formula.clone(),
                    message: error.to_string(),
                }
            })?;
            if parsed.is_empty() {
                return Err(ElementInventoryError::MalformedFormula {
                    formula: carrier.formula,
                    message: "formula contains no elements".into(),
                });
            }
            for (element, count) in parsed {
                let element = normalize_element(element)?;
                let contribution = carrier.amount * count as f64;
                validate_amount(&element, contribution)?;
                *amounts.entry(element).or_insert(0.0) += contribution;
            }
        }
        Self::from_canonical_amounts(amounts)
    }

    /// Returns the amount for one canonical element symbol, if present.
    pub fn amount(&self, element: &str) -> Option<f64> {
        self.amounts.get(element).copied()
    }

    /// Returns canonical element/amount pairs in stable order.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&str, f64)> {
        self.amounts
            .iter()
            .map(|(element, amount)| (element.as_str(), *amount))
    }

    /// Returns canonical element labels in the same order as [`Self::entries`].
    pub fn element_labels(&self) -> impl ExactSizeIterator<Item = &str> {
        self.amounts.keys().map(String::as_str)
    }

    pub fn len(&self) -> usize {
        self.amounts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.amounts.is_empty()
    }

    /// Aligns `b` to an already canonical matrix-column order.
    ///
    /// This is intentionally a pure mapping. A caller that needs a full
    /// representability check must perform it after the real species universe
    /// and its matrix `A` have been resolved.
    pub fn aligned_to(&self, element_labels: &[String]) -> Result<Vec<f64>, ElementInventoryError> {
        let known: std::collections::BTreeSet<_> =
            element_labels.iter().map(String::as_str).collect();
        if let Some(element) = self
            .amounts
            .keys()
            .find(|element| !known.contains(element.as_str()))
        {
            return Err(ElementInventoryError::ElementMissingFromLayout {
                element: element.clone(),
            });
        }
        Ok(element_labels
            .iter()
            .map(|element| self.amount(element).unwrap_or(0.0))
            .collect())
    }

    fn from_canonical_amounts(
        mut amounts: BTreeMap<String, f64>,
    ) -> Result<Self, ElementInventoryError> {
        amounts.retain(|_, amount| *amount != 0.0);
        if amounts.is_empty() {
            return Err(ElementInventoryError::EmptyInventory);
        }
        if let Some((element, amount)) = amounts
            .iter()
            .find(|(_, amount)| !amount.is_finite() || **amount <= 0.0)
        {
            return Err(ElementInventoryError::InvalidAmount {
                element: element.clone(),
                amount: *amount,
            });
        }
        Ok(Self { amounts })
    }
}

fn normalize_element(element: String) -> Result<String, ElementInventoryError> {
    let element = element.trim().to_string();
    if element.is_empty() || !is_known_element(&element) {
        return Err(ElementInventoryError::UnknownElement { element });
    }
    Ok(element)
}

fn validate_amount(subject: &str, amount: f64) -> Result<(), ElementInventoryError> {
    if !amount.is_finite() || amount < 0.0 {
        return Err(ElementInventoryError::InvalidAmount {
            element: subject.into(),
            amount,
        });
    }
    Ok(())
}

/// Rejections at the public elemental-inventory boundary.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ElementInventoryError {
    #[error("element inventory contains an unknown element '{element}'")]
    UnknownElement { element: String },
    #[error(
        "element inventory amount for '{element}' must be finite and non-negative, got {amount}"
    )]
    InvalidAmount { element: String, amount: f64 },
    #[error("formal elemental carrier formula must not be empty")]
    EmptyFormula,
    #[error("cannot parse formal elemental carrier '{formula}': {message}")]
    MalformedFormula { formula: String, message: String },
    #[error("element inventory must contain at least one positive amount")]
    EmptyInventory,
    #[error("element '{element}' is absent from the resolved matrix layout")]
    ElementMissingFromLayout { element: String },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formal_carriers_build_only_a_canonical_element_inventory() {
        let inventory = ElementInventory::from_formal_carriers([
            FormalElementCarrier::new("C2H4", 2.0).unwrap(),
            FormalElementCarrier::new("O2", 3.0).unwrap(),
        ])
        .unwrap();

        assert_eq!(
            inventory.entries().collect::<Vec<_>>(),
            vec![("C", 4.0), ("H", 8.0), ("O", 6.0)]
        );
        assert!(inventory.amount("C2H4").is_none());
    }

    #[test]
    fn duplicate_formal_carriers_reject_overflow_after_aggregation() {
        let error = ElementInventory::from_formal_carriers([
            FormalElementCarrier::new("H", f64::MAX).unwrap(),
            FormalElementCarrier::new("H", f64::MAX).unwrap(),
        ])
        .expect_err("the aggregate elemental amount must not become infinite");

        assert!(matches!(
            error,
            ElementInventoryError::InvalidAmount { element, amount }
                if element == "H" && amount.is_infinite()
        ));
    }

    #[test]
    fn direct_entries_merge_duplicates_remove_zeroes_and_ignore_order() {
        let left = ElementInventory::from_amounts([("O", 0.0), ("H", 2.0), ("O", 1.0)]).unwrap();
        let right = ElementInventory::from_amounts([("O", 1.0), ("H", 2.0)]).unwrap();
        assert_eq!(left, right);
        assert_eq!(
            left.entries().collect::<Vec<_>>(),
            vec![("H", 2.0), ("O", 1.0)]
        );
    }

    #[test]
    fn invalid_formal_and_direct_input_is_rejected_before_problem_assembly() {
        assert!(matches!(
            ElementInventory::from_amounts([("Xx", 1.0)]),
            Err(ElementInventoryError::UnknownElement { .. })
        ));
        assert!(matches!(
            ElementInventory::from_amounts([("H", -1.0)]),
            Err(ElementInventoryError::InvalidAmount { .. })
        ));
        assert!(matches!(
            ElementInventory::from_formal_carriers(
                [FormalElementCarrier::new("C2$", 1.0).unwrap()]
            ),
            Err(ElementInventoryError::MalformedFormula { .. })
        ));
        assert!(matches!(
            ElementInventory::from_amounts([("H", 0.0)]),
            Err(ElementInventoryError::EmptyInventory)
        ));
    }

    #[test]
    fn alignment_rejects_elements_missing_from_the_real_species_layout() {
        let inventory = ElementInventory::from_amounts([("C", 1.0), ("H", 4.0)]).unwrap();
        assert_eq!(
            inventory.aligned_to(&["H".into(), "C".into()]).unwrap(),
            vec![4.0, 1.0]
        );
        assert!(matches!(
            inventory.aligned_to(&["H".into()]),
            Err(ElementInventoryError::ElementMissingFromLayout { .. })
        ));
    }
}
