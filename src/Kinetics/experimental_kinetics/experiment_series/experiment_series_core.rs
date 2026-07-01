//! Ядро контейнера `TGASeries`.
//!
//! Здесь живут индекс, доступ к экспериментам по id, базовые операции
//! добавления/удаления и служебные методы, от которых зависит вся серия.

use super::*;

impl TGASeries {
    pub fn new() -> Self {
        Self {
            experiments: Vec::new(),
            exp_map: HashMap::new(),
        }
    }

    pub fn push_from_file(&mut self, path: &Path) -> Result<(), TGADomainError> {
        let exp = TGAExperiment::from_csv_universal(path)?;
        self.push(exp)?;
        Ok(())
    }

    pub fn rebuild_index(&mut self) {
        self.exp_map.clear();
        for (idx, exp) in self.experiments.iter().enumerate() {
            self.exp_map.insert(exp.meta.id.clone(), idx);
        }
    }

    fn has_experiment_id(&self, id: &str) -> bool {
        self.exp_map.contains_key(id)
    }

    fn ensure_experiment_id_is_available(
        &self,
        id: &str,
        current_index: Option<usize>,
    ) -> Result<(), TGADomainError> {
        match self.has_experiment_id(id) {
            false => Ok(()),
            true => match self.exp_map.get(id).copied() {
                Some(existing) if Some(existing) == current_index => Ok(()),
                _ => Err(TGADomainError::DuplicateExperimentId(id.to_string())),
            },
        }
    }

    pub fn drop_experiment(&mut self, idx: usize) -> Result<(), TGADomainError> {
        if idx >= self.experiments.len() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' out of range",
                idx
            )));
        }
        self.experiments.remove(idx);
        self.rebuild_index();
        Ok(())
    }

    pub fn drop_experiment_by_id(&mut self, id: &str) -> Result<(), TGADomainError> {
        let idx = self
            .exp_map
            .get(id)
            .copied()
            .ok_or_else(|| TGADomainError::ExperimentNotFound(id.to_string()))?;
        self.experiments.remove(idx);
        self.rebuild_index();
        Ok(())
    }

    pub fn push(&mut self, exp: TGAExperiment) -> Result<(), TGADomainError> {
        self.ensure_experiment_id_is_available(&exp.meta.id, None)?;
        self.experiments.push(exp);
        self.rebuild_index();
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.experiments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.experiments.is_empty()
    }

    pub fn experiments(&self) -> &[TGAExperiment] {
        &self.experiments
    }

    pub fn experiments_mut(&mut self) -> &mut [TGAExperiment] {
        &mut self.experiments
    }

    pub fn get_experiment(&self, index: usize) -> Option<&TGAExperiment> {
        self.experiments.get(index)
    }

    pub fn get_experiment_mut(&mut self, index: usize) -> Option<&mut TGAExperiment> {
        self.experiments.get_mut(index)
    }

    pub fn get_experiment_by_id(&self, id: &str) -> Result<&TGAExperiment, TGADomainError> {
        self.exp_map
            .get(id)
            .and_then(|&idx| self.experiments.get(idx))
            .ok_or_else(|| TGADomainError::ExperimentNotFound(id.to_string()))
    }

    pub fn get_experiment_by_id_mut(
        &mut self,
        id: &str,
    ) -> Result<&mut TGAExperiment, TGADomainError> {
        let idx = self.index_by_id(id)?;
        self.experiments
            .get_mut(idx)
            .ok_or_else(|| TGADomainError::ExperimentNotFound(id.to_string()))
    }

    pub fn index_by_id(&self, id: &str) -> Result<usize, TGADomainError> {
        let idx = *self
            .exp_map
            .get(id)
            .ok_or_else(|| TGADomainError::ExperimentNotFound(id.to_string()))?;
        Ok(idx)
    }

    pub fn set_experiment_id(
        &mut self,
        index: usize,
        new_id: impl Into<String>,
    ) -> Result<(), TGADomainError> {
        let new_id = new_id.into();
        if self.experiments.get(index).is_none() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Experiment index '{}' out of range",
                index
            )));
        }
        self.ensure_experiment_id_is_available(&new_id, Some(index))?;
        let exp = self.experiments.get_mut(index).ok_or_else(|| {
            TGADomainError::InvalidOperation(format!("Experiment index '{}' out of range", index))
        })?;
        exp.meta.id = new_id;
        self.rebuild_index();
        Ok(())
    }

    pub fn remove_experiment_by_id(&mut self, id: &str) -> Result<TGAExperiment, TGADomainError> {
        let idx = self.index_by_id(id)?;
        let removed = self.experiments.remove(idx);
        self.rebuild_index();
        Ok(removed)
    }

    pub fn ids(&self) -> Vec<String> {
        self.experiments
            .iter()
            .map(|exp| exp.meta.id.clone())
            .collect()
    }

    /// Human-readable provenance chain for a column in a specific experiment.
    pub fn column_provenance_text(
        &self,
        id: &str,
        col: &str,
    ) -> Result<Option<String>, TGADomainError> {
        let exp = self.get_experiment_by_id(id)?;
        Ok(exp.dataset.column_provenance_text(col))
    }

    pub fn get_list_of_complex_id(&self) -> Vec<String> {
        let mut vec_of_complex_ids = Vec::new();
        for experiment in &self.experiments {
            let vec_of_columns_for_one_exp = experiment.list_of_columns();
            let complex_id = format!(
                "{}|{}",
                experiment.meta.id,
                vec_of_columns_for_one_exp.join(",")
            );
            vec_of_complex_ids.push(complex_id);
        }
        vec_of_complex_ids
    }
}
