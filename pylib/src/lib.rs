use pyo3::prelude::*;

mod kernels;

/// Attach a submodule and make it importable.
///
/// `add_submodule` alone only sets an attribute, so `import mispy.kernels` and
/// `from mispy.kernels import ...` both fail while `mispy.kernels.f()` works.
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
fn mispy(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_submodule(m, "kernels", kernels::register)?;
    Ok(())
}
