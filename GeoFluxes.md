# GeoFluxes  
**Multiphase Flow and Transport Modeling for Subsurface Energy and Storage Systems**

---

## Executive Summary

**GeoFluxes** is a comprehensive, Python-driven computational framework for subsurface flow, transport, and coupled process modeling across energy, environmental, and storage applications. The platform integrates multiple, well-established numerical simulators—including **PFLOTRAN**, **FEHM**, **MODFLOW**, **VS2DI**, and **VS2DTI**—within a unified, transparent, and reproducible workflow.

GeoFluxes is designed to support the **full modeling lifecycle**, from conceptual model development and scenario design to simulation execution, postprocessing, parameter estimation, uncertainty analysis, and decision-ready reporting. By combining solver-appropriate physics with a consistent Python orchestration layer, the framework preserves the numerical strengths of each engine while eliminating fragmentation in setup, analysis, QA/QC, and documentation.

The platform delivers **traceable simulation packages**, operational scenarios with explicit physical and engineering constraints, and **decision-ready outputs** suitable for technical review, regulatory-style documentation, and stakeholder communication. Emphasis is placed on physical realism, transparent assumptions, reproducibility, and defensible uncertainty framing, enabling robust evaluation of **groundwater flow**, **multiphase and gas–liquid systems**, **reactive transport**, **near-field and variably saturated processes**, **microfluidic and pore-scale behavior**, and **long-term system performance** across a wide range of subsurface applications.

---

## Capabilities

GeoFluxes provides an integrated, multi-scale modeling capability spanning laboratory, field, and regional domains. Each capability leverages the most appropriate numerical engines while maintaining a consistent workflow for scenario definition, execution control, postprocessing, parameter estimation, QA/QC, and reporting.

---

## Groundwater Flow and Solute Transport

GeoFluxes supports saturated groundwater flow and solute transport modeling for environmental, energy, and infrastructure applications.

### What We Can Do
- Regional- and site-scale groundwater flow modeling  
- Hydraulic head and pressure evolution under natural and engineered stresses  
- Conservative and reactive solute transport  
- Recharge, boundary-condition, and property sensitivity analysis  

**Primary engines**
- MODFLOW and MODFLOW-family tools  

**Supporting engines**
- PFLOTRAN (density-dependent flow and coupled transport)  
- FEHM (flow with thermal coupling)  

---

## Multiphase Flow (CO₂, H₂, Gas–Water Systems)

GeoFluxes supports multiphase flow and gas–liquid interaction modeling across reservoir-scale, near-field, and transition-zone regimes.

### What We Can Do
- Gas injection scenarios with pressure buildup and plume migration  
- Buoyancy-driven flow, containment, and trapping behavior  
- Salinity- and temperature-dependent effects  
- Coupled saturated–unsaturated transitions near wells and boundaries  

**Primary engines**
- PFLOTRAN — multiphase gas–water flow and compositional effects  
- FEHM — multiphase and thermal–hydrologic systems  

**Supporting engines**
- MODFLOW-family (density-aware formulations)  
- VS2DI / VS2DTI for near-field and partially saturated zones  

---

## Wells, Operations, and Engineering Constraints

GeoFluxes represents coupled wells and operational controls under realistic engineering and regulatory constraints.

### What We Can Do
- Injection and production scheduling (CO₂, H₂, brine, water)  
- Pressure, rate, and fracture-gradient constraint enforcement  
- Well interference analysis and operational trade studies  
- Field-level performance metrics and diagnostics  

**Primary engines**
- PFLOTRAN  
- FEHM  

**Supporting workflows**
- MODFLOW well packages  
- Python-based operational logic for scheduling and constraint handling  

---

## Reactive Transport and Geochemical Processes

GeoFluxes incorporates reactive transport where geochemistry governs system evolution and long-term behavior.

### What We Can Do
- Coupled flow, transport, and geochemical reaction modeling  
- Mineral dissolution and precipitation  
- Gas partitioning and phase exchange  
- Time-dependent mass balance and reaction-path analysis  
- Permanence, storage, and transformation metrics with explicit assumptions  

**Primary engines**
- PFLOTRAN  

**Supporting workflows**
- Python-based reaction bookkeeping, diagnostics, and reporting  

---

## Variably Saturated Flow and Near-Field Processes

GeoFluxes supports vadose-zone and near-surface modeling where capillary effects and partial saturation dominate.

### What We Can Do
- Variably saturated flow and transport simulations  
- Infiltration, leakage, and near-field process representation  
- Coupling of vadose-zone and deeper saturated or multiphase systems  

**Primary engines**
- VS2DI  
- VS2DTI  

---

## Microfluidic and Pore-Scale Simulations

GeoFluxes supports microfluidic, pore-scale, and laboratory-scale simulations for process understanding, benchmarking, and scale-bridging.

### What We Can Do
- Flow and transport in microfluidic and lab-on-chip geometries  
- Geometry-controlled experiments and synthetic pore networks  
- Upscaling insights to inform continuum-scale models  

**Primary engines**
- PFLOTRAN (micro-domain configurations)  
- FEHM (small-scale domains)  

**Supporting workflows**
- Python-based geometry generation, meshing, and analysis  

---

## Parameter Estimation, Inverse Modeling, and Uncertainty Analysis

GeoFluxes provides parameter estimation and inverse modeling capabilities across all supported simulators.

### What We Can Do
- Calibration of hydraulic, transport, thermal, and reactive parameters  
- Gradient-free, ensemble-based, and sampling-driven inverse workflows  
- Sensitivity analysis, identifiability assessment, and uncertainty propagation  
- Cross-solver parameter comparison and risk framing  

**Supported engines**
- PFLOTRAN  
- FEHM  
- MODFLOW  
- VS2DI  
- VS2DTI  

**Integration layer**
- Python-based optimization, sampling, and diagnostic tools  

---

## Postprocessing, KPIs, and Reporting

GeoFluxes converts raw simulator outputs into consistent, decision-ready products.

### What We Can Do
- Automated plots for pressure, flow, plume extent, and transport metrics  
- Mass balance, storage partitioning, and diagnostic KPIs  
- Cross-solver comparison and standardized reporting  
- Outputs suitable for technical, regulatory, and stakeholder review  

---

## Reproducibility, QA/QC, and Audit Readiness

GeoFluxes is designed to produce defensible, auditable results.

### What We Can Do
- Fully reproducible run packages with complete provenance  
- Solver-agnostic QA/QC (mass balance, stability, convergence checks)  
- Versioning, change tracking, and assumptions registers  
- Transparent linkage from inputs to final deliverables  

---

## Why This Matters for Clients

GeoFluxes reduces technical risk by transforming complex subsurface problems into **transparent, testable, and repeatable modeling workflows**. By supporting groundwater flow, multiphase systems, reactive transport, parameter estimation, and pore-scale processes within a single coherent platform, GeoFluxes enables faster insight, clearer uncertainty communication, and stronger technical justification for decisions related to energy systems, environmental management, and subsurface storage.

---

## License

© 2026 Md Lal Mamud. All rights reserved.

GeoFluxes is proprietary research software. Use, reproduction, modification, or distribution of this software requires explicit prior written permission from the author. See the `LICENSE` file for full terms.

---

## Citation

If you use **GeoFluxes** in research, proposals, or technical studies, please cite appropriately. Citation details for GeoFluxes will be provided in a future release.

Users are responsible for citing the underlying numerical engines (**PFLOTRAN**, **FEHM**, **MODFLOW**, **VS2DI**, **VS2DTI**) according to their respective citation guidelines and licensing requirements.
