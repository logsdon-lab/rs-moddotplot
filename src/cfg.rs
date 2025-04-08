#[derive(Debug, Clone)]
pub struct SelfIdentConfig {
    pub window_size: usize,
    pub k: usize,
    pub id_threshold: f32,
    pub delta: f32,
    pub modimizer: usize,
}

impl Default for SelfIdentConfig {
    fn default() -> Self {
        Self {
            window_size: 5000,
            k: 21,
            id_threshold: 0.86,
            delta: 0.5,
            modimizer: 1000,
        }
    }
}
