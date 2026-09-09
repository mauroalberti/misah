use pyo3::prelude::*;

mod fields;
mod kernels;
mod mechanisms;

/// Attach a submodule and make it importable.
///
/// `add_submodule` alone only sets an attribute, so `import misah.kernels` and
/// `from misah.kernels import ...` both fail while `misah.kernels.f()` works.
/// Registering it in `sys.modules` under its dotted path is what closes that gap.
fn add_submodule(
    parent: &Bound<'_, PyModule>,
    name: &str,
    build: impl FnOnce(&Bound<'_, PyModule>) -> PyResult<()>,
) -> PyResult<()> {
    let py = parent.py();

    let child = PyModule::new(py, name)?;
    build(&child)?;
    parent.add_submodule(&child)?;

    let dotted = format!("{}.{}", parent.name()?, name);
    py.import("sys")?
        .getattr("modules")?
        .set_item(dotted, &child)?;

    Ok(())
}

#[pymodule]
fn _misah(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Three Rust files, one Python module. `fields.rs` and `mechanisms.rs`
    // were split off `kernels.rs` for length, and a second import path would
    // have made a decision about source layout into something callers have to
    // know.
    add_submodule(m, "kernels", |kernels| {
        kernels::register(kernels)?;
        fields::register(kernels)?;
        mechanisms::register(kernels)
    })?;
    Ok(())
}
