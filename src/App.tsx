/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState, useEffect, useMemo } from "react";
import { 
  AppTab, 
  SimulationParams, 
  SimulationResult, 
  UploadedDataset, 
  PotentialModelType, 
  NuclearModelType, 
  ExchangeModelType, 
  ProjectileType,
  SavedProfile,
  PresetCompound
} from "./types";
import { runScatteringSimulation, PRESET_COMPOUNDS, PERIODIC_TABLE } from "./utils/physicsSolver";
import PeriodicTable from "./components/PeriodicTable";
import CustomChart from "./components/CustomChart";
import FileUploader from "./components/FileUploader";
import WasmCompilerWidget from "./components/WasmCompilerWidget";
import TheoryGuide from "./components/TheoryGuide";

import { 
  Atom, 
  FlaskConical, 
  Settings, 
  Upload, 
  Cpu, 
  BookOpen, 
  FileText, 
  Terminal, 
  Sparkles, 
  Download, 
  Undo,
  ListFilter,
  Sun,
  Moon,
  Trash2,
  Layers,
  Plus,
  Minus,
  Info,
  Layers3,
  Flame,
  Activity,
  Eye,
  EyeOff
} from "lucide-react";

export default function App() {
  // Navigation tabs
  const [activeTab, setActiveTab] = useState<AppTab>(AppTab.SIMULATOR);

  // Theme control state (Persist in standard client-side localStorage)
  const [isDarkMode, setIsDarkMode] = useState<boolean>(() => {
    const saved = localStorage.getItem("theme");
    return saved === null ? true : saved !== "light";
  });

  // Base simulation configuration modes: "single" element vs molecular "compound"
  const [simulationMode, setSimulationMode] = useState<"single" | "compound">("single");

  // Single Z target state
  const [selectedZ, setSelectedZ] = useState<number>(79); // Default Gold (Au)

  // Compound variables
  const [compoundFormula, setCompoundFormula] = useState<string>("H2O");
  const [compoundName, setCompoundName] = useState<string>("Water");
  const [compoundAtoms, setCompoundAtoms] = useState<{ symbol: string; atomicNumber: number; stoichiometry: number }[]>([
    { symbol: "H", atomicNumber: 1, stoichiometry: 2 },
    { symbol: "O", atomicNumber: 8, stoichiometry: 1 }
  ]);

  // Command/Incident beam parameters states
  const [projectile, setProjectile] = useState<ProjectileType>("electron");
  const [energy, setEnergy] = useState<number>(200); // Default 200 eV
  const [energyUnit, setEnergyUnit] = useState<"eV" | "keV" | "MeV">("eV");

  // Advanced physical approximation parameters
  const [potentialModel, setPotentialModel] = useState<PotentialModelType>("dirac-fock");
  const [nuclearModel, setNuclearModel] = useState<NuclearModelType>("fermi");
  const [exchangeModel, setExchangeModel] = useState<ExchangeModelType>("riley-truhlar");
  
  const [absorptionModel, setAbsorptionModel] = useState<boolean>(true);
  const [correlationPolarization, setCorrelationPolarization] = useState<boolean>(true);
  const [gridPoints, setGridPoints] = useState<number>(500);

  // Active scattering simulation state
  const [simulation, setSimulation] = useState<SimulationResult | null>(null);
  const [simulationEngine, setSimulationEngine] = useState<"native_fortran" | "ts_emulated">("ts_emulated");
  const [simLoading, setSimLoading] = useState<boolean>(false);
  const [customWasmSolver, setCustomWasmSolver] = useState<any>(null);
  const [wasmBinarySize, setWasmBinarySize] = useState<number>(0);

  // Saved overlays comparisons state
  const [savedProfiles, setSavedProfiles] = useState<SavedProfile[]>([]);

  // Batch Range parameter sweeps state
  const [selectedSweepEnergies, setSelectedSweepEnergies] = useState<number[]>([100, 500, 1000]);

  // Compound Builder active inputs
  const [customAtomStoichiometry, setCustomAtomStoichiometry] = useState<number>(1);

  // Custom uploaded spreadsheet experimental datasets
  const [uploadedDatasets, setUploadedDatasets] = useState<UploadedDataset[]>([]);

  // Interface view toggles
  const [showPeriodicTable, setShowPeriodicTable] = useState<boolean>(true);
  const [activeFileView, setActiveFileView] = useState<"input" | "output">("output");

  // Theme Sync effect (Applies dark class to document body elegantly)
  useEffect(() => {
    if (isDarkMode) {
      document.documentElement.classList.add("dark");
      localStorage.setItem("theme", "dark");
    } else {
      document.documentElement.classList.remove("dark");
      localStorage.setItem("theme", "light");
    }
  }, [isDarkMode]);

  // Re-run simulation instantly whenever any parameter updates
  useEffect(() => {
    let active = true;
    const params: SimulationParams = {
      mode: simulationMode,
      atomicNumber: selectedZ,
      projectile,
      energy,
      energyUnit,
      potentialModel,
      nuclearModel,
      exchangeModel,
      absorptionModel,
      correlationPolarization,
      gridPoints,
      compoundFormula,
      compoundName,
      compoundAtoms: simulationMode === "compound" ? compoundAtoms.map(a => ({ atomicNumber: a.atomicNumber, stoichiometry: a.stoichiometry })) : undefined
    };

    const triggerSimulation = async () => {
      setSimLoading(true);

      if (customWasmSolver) {
        try {
          if (simulationMode === "single") {
            const mappedInputs = {
              atomicNumber: selectedZ,
              projectile: projectile === "electron" ? 1 : -1,
              energy: energyUnit === "keV" ? energy * 1000 : energyUnit === "MeV" ? energy * 1000000 : energy,
              potentialModel: potentialModel === "dirac-fock" ? 1 : potentialModel === "hartree-fock" ? 2 : potentialModel === "bohr-screening" ? 3 : 4,
              nuclearModel: nuclearModel === "point" ? 1 : nuclearModel === "uniform" ? 2 : 3,
              exchangeModel: exchangeModel === "none" ? 0 : exchangeModel === "furness-mccarthy" ? 1 : 2,
              absorptionModel: absorptionModel ? 1 : 0,
              correlationPolarization: correlationPolarization ? 1 : 0
            };
            const wasmResult = customWasmSolver(mappedInputs);
            if (wasmResult.success && active) {
              const formattedResult = {
                params,
                element: PERIODIC_TABLE.find(el => el.number === selectedZ) || PERIODIC_TABLE[5],
                dcsData: wasmResult.dcsData,
                phaseShifts: wasmResult.phaseShifts,
                totalElasticCrossSection: wasmResult.totalElasticCrossSection,
                momentumTransferCrossSection: wasmResult.totalElasticCrossSection * 0.9,
                deBroglieWavelength: wasmResult.deBroglieWavelength,
                elsepaInFile: `#################################################################\n# ELSEPA Input file compiled via Web-Assembly (WASM Core) #\n#################################################################\n`,
                elsepaOutFile: wasmResult.elsepaOut,
                timestamp: new Date().toLocaleTimeString()
              };
              setSimulation(formattedResult);
              setSimulationEngine("native_fortran");
              setSimLoading(false);
              return;
            }
          } else {
            // Compound mode utilizing WASM solvers
            const subResults = compoundAtoms.map(atom => {
              const elData = PERIODIC_TABLE.find(el => el.number === atom.atomicNumber) || PERIODIC_TABLE[5];
              const mappedInputs = {
                atomicNumber: atom.atomicNumber,
                projectile: projectile === "electron" ? 1 : -1,
                energy: energyUnit === "keV" ? energy * 1000 : energyUnit === "MeV" ? energy * 1000000 : energy,
                potentialModel: potentialModel === "dirac-fock" ? 1 : potentialModel === "hartree-fock" ? 2 : potentialModel === "bohr-screening" ? 3 : 4,
                nuclearModel: nuclearModel === "point" ? 1 : nuclearModel === "uniform" ? 2 : 3,
                exchangeModel: exchangeModel === "none" ? 0 : exchangeModel === "furness-mccarthy" ? 1 : 2,
                absorptionModel: absorptionModel ? 1 : 0,
                correlationPolarization: correlationPolarization ? 1 : 0
              };
              const wasmResult = customWasmSolver(mappedInputs);
              return {
                stoichiometry: atom.stoichiometry,
                element: elData,
                result: wasmResult
              };
            });

            const combinedDcsData: any[] = [];
            let combinedTotalElasticXC = 0.0;
            const deBroglieWavelength = subResults[0]?.result.deBroglieWavelength || 1.0;

            for (let angleDeg = 0; angleDeg <= 180; angleDeg++) {
              let weightedDcsSum = 0.0;
              let weightedRutherfordSum = 0.0;
              let weightedShermanNumer = 0.0;

              subResults.forEach(sub => {
                const pt = sub.result.dcsData[angleDeg];
                if (pt) {
                  const wDcs = sub.stoichiometry * pt.dcs;
                  weightedDcsSum += wDcs;
                  weightedRutherfordSum += sub.stoichiometry * pt.dcsRutherford;
                  weightedShermanNumer += wDcs * pt.Sherman;
                }
              });

              combinedDcsData.push({
                angle: angleDeg,
                dcs: weightedDcsSum,
                dcsRutherford: weightedRutherfordSum,
                Sherman: weightedDcsSum > 0 ? weightedShermanNumer / weightedDcsSum : 0
              });
            }

            subResults.forEach(sub => {
              combinedTotalElasticXC += sub.stoichiometry * sub.result.totalElasticCrossSection;
            });

            const compoundElement = {
              number: Math.round(compoundAtoms.reduce((acc, a) => acc + a.atomicNumber * a.stoichiometry, 0)),
              symbol: compoundFormula || "Compound",
              name: compoundName || "Chemical Compound",
              mass: compoundAtoms.reduce((acc, a) => {
                const el = PERIODIC_TABLE.find(e => e.number === a.atomicNumber);
                return acc + (el ? el.mass : 0) * a.stoichiometry;
              }, 0),
              category: "compound",
              group: 0,
              period: 0,
              block: ""
            };

            const elsepaOutFile = `*****************************************************************\n*                                                               *\n*   ELSEPA -- MOLECULAR COMPOSITE ANALYZER (WASM CORE)          *\n*   Independent Atom Approximation elastic scattering solver    *\n*                                                               *\n*****************************************************************\n Target compound ............. ${compoundName} (${compoundFormula})\n Total virtual charges ........ Z_eff = ${compoundElement.number}\n Projectile .................. ${projectile === "electron" ? "Electrons" : "Positrons"}\n Energy ...................... ${energy} ${energyUnit}\n\n COMPOSITION DETAILS:\n${subResults.map(sub => `  * ${sub.element.symbol} (Z=${sub.element.number}) x ${sub.stoichiometry}\n    WASM sigma_el = ${sub.result.totalElasticCrossSection.toExponential(4)} a0^2`).join("\n")}\n\n INTEGRATED COMPOUND CROSS SECTIONS (WASM):\n  Total elastic cross section (sigma_el) = ${combinedTotalElasticXC.toExponential(6)} a0^2\n*****************************************************************\n`;

            if (active) {
              setSimulation({
                params,
                element: compoundElement,
                dcsData: combinedDcsData,
                phaseShifts: [],
                totalElasticCrossSection: combinedTotalElasticXC,
                momentumTransferCrossSection: combinedTotalElasticXC * 0.9,
                deBroglieWavelength,
                elsepaInFile: ``,
                elsepaOutFile,
                timestamp: new Date().toLocaleTimeString()
              });
              setSimulationEngine("native_fortran");
              setSimLoading(false);
              return;
            }
          }
        } catch (wasmErr) {
          console.error("Custom WASM execution failed:", wasmErr);
        }
      }

      try {
        const response = await fetch("/api/simulate", {
          method: "POST",
          headers: {
            "Content-Type": "application/json"
          },
          body: JSON.stringify(params)
        });

        if (!response.ok) {
          throw new Error("HTTP simulation request failed");
        }

        const data = await response.json();
        if (active) {
          setSimulation(data);
          setSimulationEngine(data.engine === "composite_iaa_emulated" ? "ts_emulated" : (data.engine || "ts_emulated"));
        }
      } catch (err) {
        console.warn("Backend /api/simulate unavailable, falling back to local JS/TS mathematical emulator.", err);
        if (active) {
          const result = runScatteringSimulation(params);
          setSimulation(result);
          setSimulationEngine("ts_emulated");
        }
      } finally {
        if (active) {
          setSimLoading(false);
        }
      }
    };

    triggerSimulation();

    return () => {
      active = false;
    };
  }, [
    simulationMode,
    selectedZ,
    projectile,
    energy,
    energyUnit,
    potentialModel,
    nuclearModel,
    exchangeModel,
    absorptionModel,
    correlationPolarization,
    gridPoints,
    compoundFormula,
    compoundName,
    compoundAtoms
  ]);

  // Handle uploaded datasets
  const handleAddDataset = (ds: UploadedDataset) => {
    setUploadedDatasets((prev) => [...prev, ds]);
  };

  const handleRemoveDataset = (id: string) => {
    setUploadedDatasets((prev) => prev.filter(ds => ds.id !== id));
  };

  const handleToggleDatasetVisibility = (id: string) => {
    setUploadedDatasets((prev) => 
      prev.map(ds => ds.id === id ? { ...ds, visible: !ds.visible } : ds)
    );
  };

  // Saved Comparison Profiles handlers
  const handleToggleProfileVisibility = (id: string) => {
    setSavedProfiles((prev) =>
      prev.map((p) => (p.id === id ? { ...p, visible: !p.visible } : p))
    );
  };

  const handleRemoveProfile = (id: string) => {
    setSavedProfiles((prev) => prev.filter((p) => p.id !== id));
  };

  const handleClearProfiles = () => {
    setSavedProfiles([]);
  };

  const handleSaveCurrentProfile = () => {
    if (!simulation) return;

    const overlayColors = ["#ef4444", "#3b82f6", "#10b981", "#f59e0b", "#ec4899", "#8b5cf6", "#06b6d4"];
    const randomColor = overlayColors[savedProfiles.length % overlayColors.length];

    const label = simulationMode === "compound"
      ? `${compoundFormula} Compound (${compoundName}) @ ${energy} ${energyUnit}`
      : `${simulation.element.symbol} (Z=${simulation.element.number}) @ ${energy} ${energyUnit}`;

    const newProfile: SavedProfile = {
      id: `profile-${Date.now()}`,
      name: label,
      color: randomColor,
      visible: true,
      params: simulation.params,
      result: simulation
    };
    setSavedProfiles((prev) => [...prev, newProfile]);
  };

  // Parameter Sweeps Workstation handlers
  const handleToggleSweepEnergy = (eSize: number) => {
    setSelectedSweepEnergies((prev) =>
      prev.includes(eSize) ? prev.filter((s) => s !== eSize) : [...prev, eSize]
    );
  };

  const handleRunEnergySweep = () => {
    setSimLoading(true);
    const resultsToSave: SavedProfile[] = [];
    const colors = ["#ef4444", "#10b981", "#3b82f6", "#f59e0b", "#ec4899", "#8b5cf6"];

    try {
      selectedSweepEnergies.forEach((targetEnergy, idx) => {
        let unit: "eV" | "keV" | "MeV" = "eV";
        let displayVal = targetEnergy;
        if (targetEnergy >= 1000000) {
          unit = "MeV";
          displayVal = targetEnergy / 1000000;
        } else if (targetEnergy >= 1000) {
          unit = "keV";
          displayVal = targetEnergy / 1000;
        }

        const params: SimulationParams = {
          mode: simulationMode,
          atomicNumber: selectedZ,
          projectile,
          energy: displayVal,
          energyUnit: unit,
          potentialModel,
          nuclearModel,
          exchangeModel,
          absorptionModel,
          correlationPolarization,
          gridPoints,
          compoundFormula,
          compoundName,
          compoundAtoms: simulationMode === "compound" ? compoundAtoms.map(a => ({ atomicNumber: a.atomicNumber, stoichiometry: a.stoichiometry })) : undefined
        };

        const simulationResult = runScatteringSimulation(params);
        resultsToSave.push({
          id: `sweep-${Date.now()}-${targetEnergy}`,
          name: simulationMode === "compound"
            ? `${compoundFormula} Compound @ ${displayVal} ${unit}`
            : `${simulationResult.element.symbol} (Z=${simulationResult.element.number}) @ ${displayVal} ${unit}`,
          color: colors[idx % colors.length],
          visible: true,
          params,
          result: simulationResult
        });
      });

      setSavedProfiles((prev) => [...prev, ...resultsToSave]);
    } catch (err) {
      console.error("Batch sweep error:", err);
    } finally {
      setSimLoading(false);
    }
  };

  // Compound builder custom handlers
  const selectedElementOnTable = useMemo(() => {
    return PERIODIC_TABLE.find(e => e.number === selectedZ) || PERIODIC_TABLE[5];
  }, [selectedZ]);

  const handleAddCustomAtom = () => {
    // Check if element is already present, if so aggregate stoichiometry
    const existingIdx = compoundAtoms.findIndex(a => a.atomicNumber === selectedElementOnTable.number);
    let updatedAtoms = [...compoundAtoms];
    if (existingIdx >= 0) {
      updatedAtoms[existingIdx].stoichiometry += customAtomStoichiometry;
    } else {
      updatedAtoms.push({
        symbol: selectedElementOnTable.symbol,
        atomicNumber: selectedElementOnTable.number,
        stoichiometry: customAtomStoichiometry
      });
    }

    setCompoundAtoms(updatedAtoms);
    recalculateFormula(updatedAtoms);
    setCustomAtomStoichiometry(1); // reset stoichiometry
  };

  const handleRemoveCompoundAtom = (atomicNum: number) => {
    const filtered = compoundAtoms.filter(a => a.atomicNumber !== atomicNum);
    setCompoundAtoms(filtered);
    recalculateFormula(filtered);
  };

  const handleSelectPreset = (preset: PresetCompound) => {
    setCompoundFormula(preset.formula);
    setCompoundName(preset.name);
    const mapped = preset.atoms.map(a => ({
      symbol: a.symbol,
      atomicNumber: a.atomicNumber,
      stoichiometry: a.stoichiometry
    }));
    setCompoundAtoms(mapped);
  };

  const recalculateFormula = (atoms: { symbol: string, stoichiometry: number }[]) => {
    const f = atoms.map(a => `${a.symbol}${a.stoichiometry > 1 ? a.stoichiometry : ""}`).join("");
    setCompoundFormula(f || "Compound");
    setCompoundName("Custom Compound Layer");
  };

  // Triggers downloading in-memory text files (elsepa.in / elsepa.out)
  const triggerTextFileDownload = (fileName: string, content: string) => {
    const blob = new Blob([content], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = fileName;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  // Reset all parameters to physical default benchmark
  const handleResetDefaults = () => {
    setSimulationMode("single");
    setSelectedZ(79);
    setProjectile("electron");
    setEnergy(200);
    setEnergyUnit("eV");
    setPotentialModel("dirac-fock");
    setNuclearModel("fermi");
    setExchangeModel("riley-truhlar");
    setAbsorptionModel(true);
    setCorrelationPolarization(true);
    setGridPoints(500);
  };

  // Dynamic Theme Style mappings for professional light/dark looks
  const sMain = isDarkMode 
    ? "bg-[#090a0f] text-[#f8fafc] dark" 
    : "bg-slate-50 text-slate-800";
    
  const sHeader = isDarkMode 
    ? "bg-[#12131a] border-[#222430]" 
    : "bg-white border-slate-200 shadow-sm text-slate-900";
    
  const sCard = isDarkMode 
    ? "bg-[#12131a] border-[#222430]" 
    : "bg-white border-slate-200 shadow-sm text-slate-800";
    
  const sInput = isDarkMode 
    ? "bg-[#090a0f] border-[#222430] text-white" 
    : "bg-slate-50 border-slate-200 text-slate-800 focus:border-[#c5a059]";
    
  const sMuted = isDarkMode 
    ? "text-[#94a3b8]" 
    : "text-slate-500 font-semibold";
    
  const sNav = isDarkMode 
    ? "bg-[#090a0f] border-[#222430]" 
    : "bg-slate-100 border-slate-200";
    
  const sNavBtnActive = isDarkMode 
    ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm" 
    : "bg-white text-[#c5a059] border border-slate-200 shadow-xs";
    
  const sNavBtnInactive = isDarkMode 
    ? "text-[#94a3b8] hover:text-white" 
    : "text-slate-600 hover:text-slate-900";

  const sTextWhite = isDarkMode
    ? "text-white"
    : "text-slate-900 font-bold";

  return (
    <div className={`min-h-screen flex flex-col font-sans antialiased transition-colors duration-200 ${sMain}`}>
      {/* Upper Navigation Header bar */}
      <header className={`${sHeader} border-b sticky top-0 z-40 px-6 py-4 flex flex-col md:flex-row md:items-center justify-between gap-4 select-none`}>
        <div className="flex items-center gap-3">
          <div className="p-2.5 bg-[#c5a059] text-[#090a0f] rounded-xl shadow-[0_0_15px_rgba(197,160,89,0.25)] flex items-center justify-center">
            <Atom className="w-6 h-6 animate-spin" style={{ animationDuration: "14s" }} />
          </div>
          <div>
            <h1 className="text-[18px] font-bold tracking-wider font-serif flex items-center gap-1.5 leading-none">
              <span className={isDarkMode ? "text-white" : "text-slate-900"}>ELSEPA</span> 
              <span className="text-[#c5a059] font-normal italic">eScatter Quantum Suite</span>
            </h1>
            <p className="text-[11px] mt-1 leading-none text-[#c5a059] font-semibold">
              Dirac Relativistic Partial-Wave Atom & Molecular Scattering Engine
            </p>
          </div>
        </div>

        {/* Dynamic Navigation Tabs and Theme HUD */}
        <div className="flex items-center gap-3 self-end md:self-auto">
          <nav className={`flex ${sNav} p-1 rounded-xl text-xs gap-0.5`}>
            <button
              id="tab-simulator"
              onClick={() => setActiveTab(AppTab.SIMULATOR)}
              className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
                activeTab === AppTab.SIMULATOR ? sNavBtnActive : sNavBtnInactive
              }`}
            >
              <FlaskConical className="w-3.5 h-3.5" />
              Simulations HUD
            </button>

            <button
              id="tab-plotter"
              onClick={() => setActiveTab(AppTab.PLOTTER)}
              className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
                activeTab === AppTab.PLOTTER ? sNavBtnActive : sNavBtnInactive
              }`}
            >
              <Upload className="w-3.5 h-3.5" />
              Spreadsheet Plotter
            </button>

            <button
              id="tab-wasm"
              onClick={() => setActiveTab(AppTab.WASM_COMPILER)}
              className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
                activeTab === AppTab.WASM_COMPILER ? sNavBtnActive : sNavBtnInactive
              }`}
            >
              <Cpu className="w-3.5 h-3.5" />
              WASM Compiler & Core
            </button>

            <button
              id="tab-guide"
              onClick={() => setActiveTab(AppTab.PHYSICS_GUIDE)}
              className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
                activeTab === AppTab.PHYSICS_GUIDE ? sNavBtnActive : sNavBtnInactive
              }`}
            >
              <BookOpen className="w-3.5 h-3.5" />
              Theory Guide
            </button>
          </nav>

          {/* Superior Light/Dark Theme selector */}
          <button
            onClick={() => setIsDarkMode(!isDarkMode)}
            className={`p-2 rounded-xl border cursor-pointer transition flex items-center justify-center ${
              isDarkMode 
                ? "bg-[#1c1d26] border-[#c2a36b]/30 text-[#c5a059] hover:bg-[#252733]" 
                : "bg-slate-100 border-slate-250 text-amber-500 hover:bg-slate-200"
            }`}
            title={isDarkMode ? "Enable clean Light Theme" : "Enable ambient Dark Theme"}
          >
            {isDarkMode ? <Sun className="w-4.5 h-4.5" /> : <Moon className="w-4.5 h-4.5" />}
          </button>
        </div>
      </header>

      {/* Main workspace platform */}
      <main className="flex-1 p-6 flex flex-col gap-6 max-w-7xl mx-auto w-full self-center">
        {/* ACTIVE TAB 1: SIMULATOR WORKSPACES */}
        {activeTab === AppTab.SIMULATOR && (
          <div className="grid grid-cols-1 xl:grid-cols-12 gap-6 items-start">
            {/* Left Parameter Panel: Accordion controls */}
            <div className={`xl:col-span-3 flex flex-col gap-5 border rounded-2xl p-5 shadow-lg select-none ${sCard}`}>
              <div className={`flex items-center justify-between border-b pb-3 ${isDarkMode ? "border-[#222430]" : "border-slate-100"}`}>
                <h3 className={`text-sm font-bold font-serif tracking-wider flex items-center gap-2 ${sTextWhite}`}>
                  <Settings className="w-4.5 h-4.5 text-[#c5a059]" />
                  Physics controls
                </h3>
                <button
                  id="reset-defaults-lnk"
                  onClick={handleResetDefaults}
                  className="text-[10px] text-[#c5a059] hover:text-[#dfba73] font-bold uppercase tracking-widest flex items-center gap-1 cursor-pointer transition-colors"
                >
                  <Undo className="w-3 h-3" />
                  Reset Defaults
                </button>
              </div>

              {/* Medium Mode selector tabbed controller */}
              <div>
                <label className={`block text-xs font-bold mb-1.5 ${sMuted}`}>
                  Target Medium Mode
                </label>
                <div className={`grid grid-cols-2 p-0.5 border rounded-lg text-xs leading-none ${isDarkMode ? "bg-[#090a0f] border-[#222430]" : "bg-slate-100 border-slate-200"}`}>
                  <button
                    onClick={() => setSimulationMode("single")}
                    className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                      simulationMode === "single" ? sNavBtnActive : sNavBtnInactive
                    }`}
                  >
                    Atomic Element
                  </button>
                  <button
                    onClick={() => setSimulationMode("compound")}
                    className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                      simulationMode === "compound" ? sNavBtnActive : sNavBtnInactive
                    }`}
                  >
                    Chemical Compound
                  </button>
                </div>
              </div>

              {/* Compound builder workspace panel */}
              {simulationMode === "compound" && (
                <div className={`border-t pt-4 flex flex-col gap-3.5 ${isDarkMode ? "border-[#222430]" : "border-slate-100"}`}>
                  <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono flex items-center justify-between">
                    <span>Compound Workspace</span>
                    <Flame className="w-3.5 h-3.5 animate-pulse text-[#c5a059]" />
                  </div>

                  {/* Preset Compounds dropdown select */}
                  <div>
                    <label className={`block text-xs font-semibold mb-1 ${sMuted}`}>
                      Load Preset Target Compound
                    </label>
                    <select
                      onChange={(e) => {
                        const preset = PRESET_COMPOUNDS.find(p => p.formula === e.target.value);
                        if (preset) handleSelectPreset(preset);
                      }}
                      className={`w-full text-xs rounded border p-1.5 font-bold outline-none cursor-pointer ${sInput}`}
                    >
                      <option value="">-- Custom builders --</option>
                      {PRESET_COMPOUNDS.map((c) => (
                        <option key={c.formula} value={c.formula}>
                          {c.name} ({c.formula})
                        </option>
                      ))}
                    </select>
                  </div>

                  {/* Core custom list of compound ingredients */}
                  <div className={`rounded-xl p-3 flex flex-col gap-2 ${isDarkMode ? "bg-[#090a0f] border border-[#222430]" : "bg-slate-50 border border-slate-150"}`}>
                    <span className="text-[10px] font-extrabold text-[#c5a059] uppercase font-mono tracking-wider block">
                      Composition ingredients:
                    </span>
                    {compoundAtoms.length === 0 ? (
                      <span className="text-xs italic text-gray-400">Empty composite context. Click table to add.</span>
                    ) : (
                      <div className="flex flex-wrap gap-1.5 max-h-24 overflow-y-auto">
                        {compoundAtoms.map((atom) => (
                          <div 
                            key={atom.atomicNumber} 
                            className="flex items-center gap-1.5 text-xs bg-[#c5a059]/20 border border-[#c5a059]/40 text-[#c5a059] dark:text-[#dfba73] px-2.5 py-1 rounded-lg"
                          >
                            <span className="font-bold">{atom.symbol}<sub>{atom.stoichiometry}</sub></span>
                            <button 
                              onClick={() => handleRemoveCompoundAtom(atom.atomicNumber)} 
                              className="text-red-500 hover:text-red-400 font-extrabold cursor-pointer ml-1 text-[11px]"
                              title="Delete element"
                            >
                              ✕
                            </button>
                          </div>
                        ))}
                      </div>
                    )}

                    {/* Quick Add Interface incorporating Selected Element on table */}
                    <div className="border-t border-dashed border-[#c5a059]/30 mt-1.5 pt-2 flex flex-col gap-2">
                      <div className="flex items-center justify-between text-[11px] font-semibold text-[#8a99ad]">
                        <span>Selected from table:</span>
                        <strong className="text-white font-mono">{selectedElementOnTable.symbol} (Z={selectedElementOnTable.number})</strong>
                      </div>
                      <div className="flex items-center gap-2">
                        <div className="flex-1 flex items-center justify-between border border-[#222430] rounded bg-[#090a0f]/50 p-1">
                          <button 
                            onClick={() => setCustomAtomStoichiometry(prev => Math.max(1, prev - 1))}
                            className="bg-slate-800 text-white w-4.5 h-4.5 flex items-center justify-center rounded cursor-pointer"
                          >
                            -
                          </button>
                          <span className="text-xs font-mono font-bold text-white">{customAtomStoichiometry}</span>
                          <button 
                            onClick={() => setCustomAtomStoichiometry(prev => prev + 1)}
                            className="bg-slate-800 text-white w-4.5 h-4.5 flex items-center justify-center rounded cursor-pointer"
                          >
                            +
                          </button>
                        </div>
                        <button
                          onClick={handleAddCustomAtom}
                          className="px-2.5 py-1.5 bg-[#c5a059] text-[#090a0f] font-bold rounded text-xs hover:bg-[#dfba73] transition cursor-pointer"
                        >
                          Add Atom
                        </button>
                      </div>
                    </div>
                  </div>

                  {/* Formula and composite name fields */}
                  <div className="grid grid-cols-2 gap-2 text-xs">
                    <div>
                      <label className={`block text-[10px] uppercase font-bold mb-0.5 ${sMuted}`}>Formula</label>
                      <input 
                        type="text" 
                        value={compoundFormula} 
                        onChange={(e) => setCompoundFormula(e.target.value)}
                        className={`w-full p-1 border rounded text-xs text-center font-mono font-bold ${sInput}`}
                      />
                    </div>
                    <div>
                      <label className={`block text-[10px] uppercase font-bold mb-0.5 ${sMuted}`}>Compound Name</label>
                      <input 
                        type="text" 
                        value={compoundName} 
                        onChange={(e) => setCompoundName(e.target.value)}
                        className={`w-full p-1 border rounded text-xs text-center font-serif font-bold ${sInput}`}
                      />
                    </div>
                  </div>
                </div>
              )}

              {/* Incident particle beam accordion details */}
              <div className="flex flex-col gap-4 border-t border-[#222430] pt-4">
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  1. Collision Projectile
                </div>

                {/* Particle selector */}
                <div>
                  <label className={`block text-xs font-semibold mb-1.5 ${sMuted}`}>
                    Projectile Type
                  </label>
                  <div className={`grid grid-cols-2 p-0.5 border rounded-lg text-xs leading-none ${isDarkMode ? "bg-[#090a0f] border-[#222430]" : "bg-slate-50 border-slate-200"}`}>
                    <button
                      id="opt-particle-electron"
                      onClick={() => setProjectile("electron")}
                      className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                        projectile === "electron" ? sNavBtnActive : sNavBtnInactive
                      }`}
                    >
                      Electron (e⁻)
                    </button>
                    <button
                      id="opt-particle-positron"
                      onClick={() => setProjectile("positron")}
                      className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                        projectile === "positron" ? sNavBtnActive : sNavBtnInactive
                      }`}
                    >
                      Positron (e⁺)
                    </button>
                  </div>
                </div>

                {/* Beam energy sliders and precise inputs */}
                <div className="flex flex-col gap-2">
                  <div className="flex items-center justify-between">
                    <label className={`block text-xs font-semibold ${sMuted}`}>
                      Kinetic Energy
                    </label>
                    <div className={`flex p-0.5 rounded-md text-[10px] font-extrabold border ${isDarkMode ? "bg-[#090a0f] border-[#222430]" : "bg-slate-100 border-slate-250"}`}>
                      {(["eV", "keV", "MeV"] as const).map(unit => (
                        <button
                          key={unit}
                          id={`unit-${unit}`}
                          onClick={() => {
                            setEnergyUnit(unit);
                            if (unit === "eV") setEnergy(200);
                            if (unit === "keV") setEnergy(50);
                            if (unit === "MeV") setEnergy(2);
                          }}
                          className={`px-1.5 py-0.5 rounded cursor-pointer transition ${
                            energyUnit === unit ? "bg-[#c5a059] text-[#090a0f] font-bold" : "text-gray-500"
                          }`}
                        >
                          {unit}
                        </button>
                      ))}
                    </div>
                  </div>

                  <div className="flex items-center gap-2">
                    <input
                      id="input-energy-slider"
                      type="range"
                      min={energyUnit === "eV" ? "10" : energyUnit === "keV" ? "1" : "0.1"}
                      max={energyUnit === "eV" ? "2000" : energyUnit === "keV" ? "500" : "15"}
                      step={energyUnit === "eV" ? "10" : energyUnit === "keV" ? "5" : "0.1"}
                      value={energy}
                      onChange={(e) => setEnergy(parseFloat(e.target.value))}
                      className="flex-1 accent-[#c5a059] cursor-pointer h-1 bg-[#222430] rounded-lg outline-none"
                    />
                    <input
                      id="input-energy"
                      type="number"
                      value={energy}
                      onChange={(e) => setEnergy(Math.max(0.01, parseFloat(e.target.value) || 0))}
                      className={`w-16 text-center text-xs border p-1 rounded font-mono font-bold outline-none ${sInput}`}
                    />
                  </div>
                </div>
              </div>

              {/* Advanced Dirac formulation elements */}
              <div className={`flex flex-col gap-4 border-t pt-4 ${isDarkMode ? "border-[#222430]" : "border-slate-100"}`}>
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  2. Quantum Dirac Potentials
                </div>

                {/* Potentials drop blocks */}
                <div>
                  <label className={`block text-xs font-semibold mb-1 ${sMuted}`}>
                    Static Potential model
                  </label>
                  <select
                    id="select-potential-model"
                    value={potentialModel}
                    onChange={(e) => setPotentialModel(e.target.value as any)}
                    className={`w-full text-xs rounded border p-1.5 font-bold outline-none cursor-pointer ${sInput}`}
                  >
                    <option value="dirac-fock">Dirac-Fock (exact wavefunctions)</option>
                    <option value="hartree-fock">Hartree-Fock (SCF approximation)</option>
                    <option value="bohr-screening">Bohr Exponential Screening</option>
                    <option value="yukawa">Yukawa Potential</option>
                  </select>
                </div>

                {/* Local exchange models for projectile-electron overlap */}
                {projectile === "electron" && (
                  <div>
                    <label className={`block text-xs font-semibold mb-1 ${sMuted}`}>
                      Exchange Interaction forces
                    </label>
                    <select
                      id="select-exchange-model"
                      value={exchangeModel}
                      onChange={(e) => setExchangeModel(e.target.value as any)}
                      className={`w-full text-xs rounded border p-1.5 font-bold outline-none cursor-pointer ${sInput}`}
                    >
                      <option value="none">No Exchange Forces</option>
                      <option value="furness-mccarthy">Furness-McCarthy (FM model)</option>
                      <option value="riley-truhlar">Riley-Truhlar (RT Local force)</option>
                    </select>
                  </div>
                )}

                {/* Finite size charge density parameters */}
                <div>
                  <label className={`block text-xs font-semibold mb-1 ${sMuted}`}>
                    Nuclear Charge Density distribution
                  </label>
                  <select
                    id="select-nuclear-model"
                    value={nuclearModel}
                    onChange={(e) => setNuclearModel(e.target.value as any)}
                    className={`w-full text-xs rounded border p-1.5 font-bold outline-none cursor-pointer ${sInput}`}
                  >
                    <option value="point">Point Charge approximation</option>
                    <option value="uniform">Uniform Charged Sphere</option>
                    <option value="fermi">Fermi Two-parameter distribution</option>
                  </select>
                </div>
              </div>

              {/* Secondary auxiliary forces checkbox elements */}
              <div className={`flex flex-col gap-3 border-t pt-4 ${isDarkMode ? "border-[#222430]" : "border-slate-100"}`}>
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  3. Dynamic Auxiliary Potentials
                </div>

                <label className="flex items-center gap-2.5 cursor-pointer text-xs font-semibold leading-none">
                  <input
                    id="check-correlation"
                    type="checkbox"
                    checked={correlationPolarization}
                    onChange={(e) => setCorrelationPolarization(e.target.checked)}
                    className="w-4 h-4 rounded text-[#c5a059] focus:ring-0 focus:ring-offset-0 accent-[#c5a059] cursor-pointer"
                  />
                  Correlation-Polarization
                </label>

                <label className="flex items-center gap-2.5 cursor-pointer text-xs font-semibold leading-none">
                  <input
                    id="check-absorption"
                    type="checkbox"
                    checked={absorptionModel}
                    onChange={(e) => setAbsorptionModel(e.target.checked)}
                    className="w-4 h-4 rounded text-[#c5a059] focus:ring-0 focus:ring-offset-0 accent-[#c5a059] cursor-pointer"
                  />
                  Inelastic Absorption Potential
                </label>
              </div>

              {/* Dynamic live simulation metadata card */}
              {simulation && (
                <div className={`rounded-xl p-3 flex flex-col gap-1.5 text-[10.5px] leading-snug border ${isDarkMode ? "bg-[#090a0f] border-[#222430]" : "bg-slate-50 border-slate-150"}`}>
                  <div className={`font-mono font-bold border-b pb-1 mb-1 text-[9px] tracking-widest flex items-center justify-between ${isDarkMode ? "border-[#222430] text-gray-400" : "border-slate-200 text-slate-500"}`}>
                    <span>SIMULATOR RUN STATISTICS:</span>
                  </div>
                  <div className={`flex justify-between items-center px-1.5 py-1 rounded-lg border mb-1.5 select-none text-[9.5px] ${isDarkMode ? "bg-[#161720] border-[#222430]" : "bg-white border-slate-150 shadow-xs"}`}>
                    <span className="font-bold text-slate-400 font-mono">SOLVER STATUS:</span>
                    {simLoading ? (
                      <span className="text-yellow-500 font-extrabold font-mono animate-pulse uppercase flex items-center gap-1">
                        Computing...
                      </span>
                    ) : simulationEngine === "native_fortran" ? (
                      <span className="text-emerald-500 font-extrabold font-mono uppercase flex items-center gap-1">
                        Fortran 90 Core
                      </span>
                    ) : (
                      <span className="text-amber-500 font-extrabold font-mono uppercase">
                        JS/TS Emulator
                      </span>
                    )}
                  </div>
                  <div className="flex justify-between">
                    <span className="text-slate-400">de Broglie Wavelength:</span>
                    <strong className="text-[#c5a059] font-mono">{simulation.deBroglieWavelength.toFixed(5)} Å</strong>
                  </div>
                  <div className="flex justify-between items-center">
                    <span className="text-slate-400">Total Elastic σ<sub>el</sub>:</span>
                    <strong className="font-mono truncate max-w-[130px]" title={`${simulation.totalElasticCrossSection.toExponential(4)} a0²`}>
                      {simulation.totalElasticCrossSection.toExponential(4)} a₀²
                    </strong>
                  </div>
                  <div className="flex justify-between items-center">
                    <span className="text-slate-400">Momentum Transfer σ<sub>tr</sub>:</span>
                    <strong className="font-mono truncate max-w-[130px]" title={`${simulation.momentumTransferCrossSection.toExponential(4)} a0²`}>
                      {simulation.momentumTransferCrossSection.toExponential(4)} a₀²
                    </strong>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-slate-400">Partial Waves (Lmax):</span>
                    <strong className="font-mono">{simulation.phaseShifts.length > 0 ? simulation.phaseShifts.length - 1 : 120}</strong>
                  </div>
                  <div className="border-t border-dashed mt-1.5 pt-1.5 flex flex-col gap-1.5 border-[#383a4d]">
                    <button
                      onClick={handleSaveCurrentProfile}
                      className="w-full text-center py-1 bg-amber-500/10 border border-amber-500/30 text-amber-500 dark:text-amber-400 hover:bg-amber-500 font-bold rounded-lg hover:text-[#090a0f] transition h-7 text-[10.5px] cursor-pointer"
                    >
                      Overlay Active Curve
                    </button>
                  </div>
                </div>
              )}
            </div>

            {/* Center + Right workspace modules */}
            <div className="xl:col-span-9 flex flex-col gap-6">
              {/* Element periodic mapper elements (hide/show trigger drawer) */}
              {simulationMode === "single" && (
                <div className="flex flex-col gap-2">
                  <div className={`flex items-center justify-between border rounded-xl px-4 py-3 shadow-sm select-none ${sCard}`}>
                    <span className="text-xs font-bold flex items-center gap-1.5 font-mono leading-none">
                      <ListFilter className="w-4 h-4 text-[#c5a059] shrink-0 font-bold" />
                      Interactive Periodic Element Mapper
                    </span>
                    <button
                      id="table-toggle-btn"
                      onClick={() => setShowPeriodicTable(!showPeriodicTable)}
                      className="text-xs text-[#c5a059] font-extrabold hover:text-[#dfba73] cursor-pointer transition-colors"
                    >
                      {showPeriodicTable ? "Collapse periodic grid" : "Expand periodic grid"}
                    </button>
                  </div>
                  {showPeriodicTable && (
                    <PeriodicTable
                      selectedZ={selectedZ}
                      onSelectElement={(atomicNumber) => setSelectedZ(atomicNumber)}
                    />
                  )}
                </div>
              )}

              {/* Dynamic Parameter sweeping module widget */}
              <div className={`border rounded-2xl p-5 shadow-md flex flex-col gap-4 ${sCard}`}>
                <div className="flex items-center gap-2">
                  <Layers className="w-4.5 h-4.5 text-[#c5a059]" />
                  <h4 className="font-serif text-sm font-bold tracking-wider">Multi-Parametric Energy Sweeps</h4>
                </div>
                <p className="text-xs text-gray-400">
                  Enable rapid energy sweep calculation parameters to overlay and analyze multiple energy outputs simultaneously:
                </p>
                <div className="flex flex-wrap items-center gap-4 border border-[#c5a059]/10 p-3 rounded-lg bg-[rgba(197,160,89,0.02)]">
                  <div className="flex flex-wrap gap-2.5">
                    {[50, 100, 200, 500, 1000, 2000, 5000].map((energyVal) => (
                      <label 
                        key={energyVal} 
                        className={`px-2 py-1 text-xs font-mono font-bold rounded-md border flex items-center gap-1.5 cursor-pointer transition ${
                          selectedSweepEnergies.includes(energyVal)
                            ? "bg-[#c5a059]/20 border-[#c5a059] text-amber-500"
                            : "bg-[#090a0f]/40 border-slate-700/50 text-[#88909c]"
                        }`}
                      >
                        <input
                          type="checkbox"
                          checked={selectedSweepEnergies.includes(energyVal)}
                          onChange={() => handleToggleSweepEnergy(energyVal)}
                          className="w-3 h-3 text-[#c5a059] focus:ring-0 outline-none accent-[#c5a059]"
                        />
                        {energyVal >= 1000 ? `${energyVal / 1000} keV` : `${energyVal} eV`}
                      </label>
                    ))}
                  </div>
                  <button
                    onClick={handleRunEnergySweep}
                    className="ml-auto px-4 py-1.5 text-xs bg-[#c5a059] text-[#090a0f] font-extrabold rounded-lg hover:bg-[#dfba73] hover:scale-101 transition flex items-center gap-1 max-w-[190px] shadow-sm cursor-pointer whitespace-nowrap"
                  >
                    <Activity className="w-3.5 h-3.5" />
                    Compute Sweep Overlay
                  </button>
                </div>
              </div>

              {/* High-Performance interactive scientific charts plotter */}
              <CustomChart
                dcsData={simulation?.dcsData || []}
                elementSymbol={simulationMode === "compound" ? compoundFormula : (simulation?.element.symbol || "Au")}
                energyText={`${energy} ${energyUnit}`}
                uploadedDatasets={uploadedDatasets}
                onToggleDatasetVisibility={handleToggleDatasetVisibility}
                savedProfiles={savedProfiles}
                onToggleProfileVisibility={handleToggleProfileVisibility}
                onRemoveProfile={handleRemoveProfile}
              />

              {/* Live active overlay comparison tables */}
              {savedProfiles.length > 0 && (
                <div className={`border rounded-2xl p-5 shadow-md flex flex-col gap-4 ${sCard}`}>
                  <div className="flex items-center justify-between">
                    <div className="flex items-center gap-2">
                      <Layers3 className="w-4.5 h-4.5 text-[#c5a059]" />
                      <h4 className="font-serif text-sm font-bold tracking-wider">Overlaid Comparisons Ledger</h4>
                    </div>
                    <button
                      onClick={handleClearProfiles}
                      className="px-2.5 py-1 text-xs bg-red-500/10 border border-red-500/20 text-red-500 rounded-md hover:bg-red-500 hover:text-white font-semibold transition cursor-pointer"
                    >
                      Clear All Comparison Layers
                    </button>
                  </div>
                  <div className="overflow-x-auto border border-slate-700/50 rounded-xl bg-black/10">
                    <table className="w-full text-xs text-left leading-normal">
                      <thead>
                        <tr className="bg-[#090a0f]/60 text-slate-400 font-mono border-b border-slate-700/50">
                          <th className="p-3">Target Layer</th>
                          <th className="p-3">Wavelength (λ)</th>
                          <th className="p-3">Total σ<sub>el</sub></th>
                          <th className="p-3">Cross transfer σ<sub>tr</sub></th>
                          <th className="p-3 text-center">Vis.</th>
                          <th className="p-3 text-right">Delete</th>
                        </tr>
                      </thead>
                      <tbody>
                        {savedProfiles.map((p) => (
                          <tr key={p.id} className="border-b border-slate-800/40 hover:bg-[#161720]/40">
                            <td className="p-3 font-serif font-bold flex items-center gap-2">
                              <span className="w-3 h-3 rounded-full shrink-0" style={{ backgroundColor: p.color }} />
                              <span style={{ color: p.color }}>{p.name}</span>
                            </td>
                            <td className="p-3 font-mono text-amber-500">{p.result.deBroglieWavelength.toFixed(5)} Å</td>
                            <td className="p-3 font-mono text-gray-300">{p.result.totalElasticCrossSection.toExponential(4)} a₀²</td>
                            <td className="p-3 font-mono text-gray-300">{p.result.momentumTransferCrossSection.toExponential(4)} a₀²</td>
                            <td className="p-3 text-center">
                              <button
                                onClick={() => handleToggleProfileVisibility(p.id)}
                                className="p-1 cursor-pointer hover:scale-105 transition hover:text-white text-gray-400"
                              >
                                {p.visible ? <Eye className="w-4 h-4 text-emerald-500" /> : <EyeOff className="w-4 h-4 text-gray-500" />}
                              </button>
                            </td>
                            <td className="p-3 text-right">
                              <button
                                onClick={() => handleRemoveProfile(p.id)}
                                className="text-red-500 hover:text-red-400 p-1 cursor-pointer"
                              >
                                <Trash2 className="w-4 h-4" />
                              </button>
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              {/* Dynamic simulated log outputs block */}
              {simulation && (
                <div className={`border rounded-2xl overflow-hidden flex flex-col h-[340px] shadow-xl ${sCard}`}>
                  <div className={`border-b px-4 py-3 flex items-center justify-between select-none shrink-0 ${isDarkMode ? "bg-[#090a0f] border-[#222430]" : "bg-slate-50 border-slate-100"}`}>
                    <div className="flex items-center gap-3">
                      <Terminal className="w-4.5 h-4.5 text-[#c5a059] shrink-0" />
                      <div className={`font-bold text-xs leading-none flex border p-0.5 rounded-lg ${isDarkMode ? "bg-[#12131a] border-[#222430]" : "bg-white border-slate-200 shadow-xs"}`}>
                        <button
                          id="btn-log-output"
                          onClick={() => setActiveFileView("output")}
                          className={`px-3 py-1 rounded-md transition cursor-pointer text-[11px] ${
                            activeFileView === "output" ? sNavBtnActive : sNavBtnInactive
                          }`}
                        >
                          elsepa.out (Tabulated Output dataset Logs)
                        </button>
                        <button
                          id="btn-log-input"
                          onClick={() => setActiveFileView("input")}
                          className={`px-3 py-1 rounded-md transition cursor-pointer text-[11px] ${
                            activeFileView === "input" ? sNavBtnActive : sNavBtnInactive
                          }`}
                        >
                          elsepa.in (Input configuration directives)
                        </button>
                      </div>
                    </div>

                    <button
                      id="download-log-file-btn"
                      onClick={() => 
                        activeFileView === "output" 
                           ? triggerTextFileDownload(`elsepa_${simulationMode === "compound" ? compoundFormula : simulation.element.symbol}_${energy}${energyUnit}.out`, simulation.elsepaOutFile)
                           : triggerTextFileDownload("elsepa.in", simulation.elsepaInFile)
                      }
                      className="text-xs bg-[#1c1d26] text-white hover:bg-[#222430] border border-[#222430] rounded-lg px-3 py-1.5 font-bold cursor-pointer transition flex items-center gap-1.5 shadow-sm"
                    >
                      <Download className="w-3.5 h-3.5 text-[#c5a059]" />
                      Download Configuration Files
                    </button>
                  </div>

                  <div className="flex-1 p-4 bg-[#040508] font-mono text-[11px] text-[#a5b4fc] overflow-auto whitespace-pre leading-normal border-t border-black select-text">
                    {activeFileView === "output" ? simulation.elsepaOutFile : simulation.elsepaInFile}
                  </div>
                </div>
              )}
            </div>
          </div>
        )}

        {/* ACTIVE TAB 2: EXCEL/CSV FILE PLOTTER DASHBOARD */}
        {activeTab === AppTab.PLOTTER && (
          <div className="grid grid-cols-1 md:grid-cols-12 gap-6 items-start">
            <div className="md:col-span-4">
              <FileUploader
                onAddDataset={handleAddDataset}
                uploadedDatasets={uploadedDatasets}
                onRemoveDataset={handleRemoveDataset}
              />
            </div>

            <div className="md:col-span-8 flex flex-col gap-5">
              <CustomChart
                dcsData={simulation?.dcsData || []}
                elementSymbol={simulationMode === "compound" ? compoundFormula : (simulation?.element.symbol || "Au")}
                energyText={`${energy} ${energyUnit}`}
                uploadedDatasets={uploadedDatasets}
                onToggleDatasetVisibility={handleToggleDatasetVisibility}
                savedProfiles={savedProfiles}
                onToggleProfileVisibility={handleToggleProfileVisibility}
                onRemoveProfile={handleRemoveProfile}
              />
            </div>
          </div>
        )}

        {/* ACTIVE TAB 3: WEB-ASSEMBLY FORTRAN COMPILER & IDE */}
        {activeTab === AppTab.WASM_COMPILER && (
          <WasmCompilerWidget 
            isDarkMode={isDarkMode}
            onCompileSuccess={(solverFn, size) => {
              setCustomWasmSolver(() => solverFn);
              setWasmBinarySize(size);
            }}
            currentZ={selectedZ}
            currentEnergy={energy}
            currentEnergyUnit={energyUnit}
          />
        )}

        {/* ACTIVE TAB 4: THEORETICAL GUIDE */}
        {activeTab === AppTab.PHYSICS_GUIDE && (
          <TheoryGuide />
        )}
      </main>

      {/* Science laboratory corporate footer info */}
      <footer className={`border-t py-6 text-center text-xs select-none font-sans font-medium mt-auto shrink-0 shadow-xs ${isDarkMode ? "bg-[#12131a] border-[#222430] text-[#94a3b8]" : "bg-white border-slate-200 text-slate-500"}`}>
        <p>© 2026 eScatter ELSEPA Computation Suite. All scattering approximations solved under relativistic quantum limits using numerical partial wave grids.</p>
      </footer>
    </div>
  );
}
