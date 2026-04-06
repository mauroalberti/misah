
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GxsurfError {
    #[error("sqlite error: {0}")]
    Sqlite(#[from] rusqlite::Error),

    #[error("missing table: {0}")]
    MissingTable(String),

    #[error("missing column '{column}' in table '{table}'")]
    MissingColumn { table: String, column: String, },

    #[error("invalid value: {0}")]
    InvalidValue(String),

    #[error("schema mismatch: {0}")]
    SchemaMismatch(String),
}
