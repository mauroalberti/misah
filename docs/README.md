# Notebooks

Worked examples of the Python bindings, in `notebooks/`. They are the long-form
counterpart to the API docstrings: the docstrings say what an argument means,
these say why you would choose one value over another and how the answer
misleads you when you choose badly.

- **`01_kernel_density.ipynb`** — kernel density estimation. What "observations
  per unit volume" buys you, why `origin` is a node and not a cell corner,
  choosing a bandwidth, why depth gets its own, and what a truncated kernel
  throws away.
- **`02_stress_field.ipynb`** — a stress tensor at every node. Why a single
  inversion over two tectonic phases returns a confident wrong answer, sizing a
  run before starting it, and the two columns that say whether a node's tensor
  is worth reading.
- **`03_focal_mechanisms.ipynb`** — the Kagan angle. Why P and T rather than
  nodal planes, why the rotation between two mechanisms has four answers and
  not one, the 120-degree bound, and turning a catalogue into a distance matrix
  to cluster it. Checked against Kagan (1991)'s own published table.

None of them reads a file. Each generates its own data, and each has something
independent to be right against: the first two build their faults through the
forward model, so what the inversion should return is known without trusting
the inversion, and the third checks itself against four numbers printed in a
1991 paper.

## Running them

```sh
pip install --pre misah        # or build the wheel: see the main README
pip install jupyter matplotlib scipy
jupyter lab docs/notebooks/
```

`matplotlib` and `scipy` are notebook dependencies, not `misah` ones. Only
`03_focal_mechanisms.ipynb` needs scipy, and only for the clustering at the
end.

## Outputs are committed

Deliberately, so the notebooks read as documentation on the web without being
run. The cost is that a diff touching a notebook is mostly base64, and the
discipline that goes with it is that outputs must not go stale: they are
regenerated whenever the code they call changes, by executing the notebook
rather than by editing the JSON.

```sh
jupyter nbconvert --to notebook --execute --inplace docs/notebooks/*.ipynb
```

The second one takes about a minute; most of it is one `stress_field` call,
which is the honest cost of the thing being demonstrated. The other two are
seconds.

## What they are not

They are not tests, and CI does not run them. The suites under `pylib/tests/`
are what assert the bindings behave — `test_fields.py` and
`test_mechanisms.py` between them cover everything these notebooks call, and do
it in a form that fails loudly. A notebook that runs is not evidence that its
numbers are right, only that nothing raised.

The prose here does make checkable claims, though, and they were checked
against what the cells actually printed rather than written from expectation.
Where the first draft and the output disagreed, the output won: the margin
discussion in the first notebook, the domain-boundary dataset in the second,
and the note in the third about why one population appears as two blobs on a
stereonet are all there because the original text claimed something the figures
did not show.
