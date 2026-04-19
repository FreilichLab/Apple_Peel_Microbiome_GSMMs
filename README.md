# Data and Code Repository: Apple Peel Microbiome Modeling

This repository serves as the data deposit for the manuscript "From Genomes to Ecological Predictions: A Standardized GSMM Curation Pipeline Uncovers *Penicillium expansum* Suppression by Native Apple Fruit Epiphytic Bacteria", [doi].

It contains the raw and curated assets used to model the metabolic interactions between the native bacterial community of the apple peel and the fungal pathogen *Penicillium expansum*. The repository is organized to provide open access to our Metagenome-Assembled Genomes (MAGs), the resulting 24 bacterial Genome-Scale Metabolic Models (GSMMs), the fungal pathogen model, and an example of custom multi-stage curation code used to generate them. Raw metagenome sequences and metabolomic data from which these models were curated and constrained are deposited in ([http://www.ncbi.nlm.nih.gov/bioproject/948808]).

---

## Repository Structure

### `MAGs/`
Contains the de-replicated sequences of Metagenome-Assembled Genomes representing the bacterial and fungal members of the apple peel community.

### `Models/`
Contains the highly curated Genome-Scale Metabolic Models (GSMMs) in standard formats (e.g., SBML/JSON).
* **Bacterial Community:** Models `GSMM_1` through `GSMM_24` and `GSMM_MR` and the curated models `[GSMM_x_]_gapfilled`
* **Fungal Pathogen:** Curated model for *Penicillium expansum*.

### `Code/`
Contains an example of the Python scripts utilized for the curation pipeline and simulation:
* **Curation Pipeline:** Scripts demonstrating the multi-stage curation workflow (diagnostic growth failure checks, iterative gap-filling, stoichiometric consistency).
* **Simulations & Analysis:** Code used to simulate community dynamics via dFBA.

---

## Contributors
* Rotem Bartuv
* Shiri Freilich
* Samir Droby

## References
Bartuv, R., Zhimo, V.Y., Sharma, V.K., Shalem, T., Droby, S., Freilich, S. 2026. "From Genomes to Ecological Predictions: A Standardized GSMM Curation Pipeline Uncovers *Penicillium expansum* Suppression by Native Apple Fruit Epiphytic Bacteria". *Microbiome* (Under review). [doi]

## Funding
This work was funded by the Israel Science Foundation [grant number 377/23].
