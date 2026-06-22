/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useState, useEffect, useMemo } from "react";
import { AppTab, SimulationParams, SimulationResult, UploadedDataset, PotentialModelType, NuclearModelType, ExchangeModelType, ProjectileType } from "./types";
import { runScatteringSimulation } from "./utils/physicsSolver";
import PeriodicTable from "./components/PeriodicTable";
import CustomChart from "./components/CustomChart";
import FileUploader from "./components/FileUploader";
import PackagerWidget from "./components/PackagerWidget";
import TheoryGuide from "./components/TheoryGuide";

import { 
  Atom, 
  FlaskConical, 
  Settings, 
  Upload, 
  Monitor, 
  BookOpen, 
  FileText, 
  Terminal, 
  HelpCircle, 
  Sparkles, 
  Download, 
  Undo,
  ListFilter
} from "lucide-react";

export default function App() {
  // Navigation tabs
  const [activeTab, setActiveTab] = useState<AppTab>(AppTab.SIMULATOR);

  // Configuration parameter states
  const [selectedZ, setSelectedZ] = useState<number>(79); // Default Gold
  const [projectile, setProjectile] = useState<ProjectileType>("electron");
  const [energy, setEnergy] = useState<number>(200); // Default 200 eV
  const [energyUnit, setEnergyUnit] = useState<"eV" | "keV" | "MeV">("eV");

  const [potentialModel, setPotentialModel] = useState<PotentialModelType>("dirac-fock");
  const [nuclearModel, setNuclearModel] = useState<NuclearModelType>("fermi");
  const [exchangeModel, setExchangeModel] = useState<ExchangeModelType>("riley-truhlar");
  
  const [absorptionModel, setAbsorptionModel] = useState<boolean>(true);
  const [correlationPolarization, setCorrelationPolarization] = useState<boolean>(true);
  const [gridPoints, setGridPoints] = useState<number>(500);

  // Computed results state
  const [simulation, setSimulation] = useState<SimulationResult | null>(null);
  const [simulationEngine, setSimulationEngine] = useState<"native_fortran" | "ts_emulated">("ts_emulated");
  const [simLoading, setSimLoading] = useState<boolean>(false);

  // Custom uploaded spreadsheet datasets
  const [uploadedDatasets, setUploadedDatasets] = useState<UploadedDataset[]>([]);

  // Open/Close sub panes
  const [showPeriodicTable, setShowPeriodicTable] = useState<boolean>(true);
  const [activeFileView, setActiveFileView] = useState<"input" | "output">("output");

  // Re-run simulation instantly whenever any primitive parameter updates
  useEffect(() => {
    let active = true;
    const params: SimulationParams = {
      atomicNumber: selectedZ,
      projectile,
      energy,
      energyUnit,
      potentialModel,
      nuclearModel,
      exchangeModel,
      absorptionModel,
      correlationPolarization,
      gridPoints
    };

    const triggerSimulation = async () => {
      setSimLoading(true);
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
          setSimulationEngine(data.engine || "ts_emulated");
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
    selectedZ,
    projectile,
    energy,
    energyUnit,
    potentialModel,
    nuclearModel,
    exchangeModel,
    absorptionModel,
    correlationPolarization,
    gridPoints
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

  return (
    <div className="min-h-screen bg-[#090a0f] text-[#f8fafc] flex flex-col font-sans antialiased text-[#e2e8f0]">
      {/* Upper Navigation Header bar */}
      <header className="bg-[#12131a] border-b border-[#222430] sticky top-0 z-40 px-6 py-4 flex flex-col md:flex-row md:items-center justify-between gap-4 select-none shadow-[0_4px_30px_rgba(0,0,0,0.4)]">
        <div className="flex items-center gap-3">
          <div className="p-2.5 bg-[#c5a059] text-[#090a0f] rounded-xl shadow-[0_0_15px_rgba(197,160,89,0.25)]">
            <Atom className="w-6 h-6 animate-spin" style={{ animationDuration: "14s" }} />
          </div>
          <div>
            <h1 className="text-[18px] font-bold tracking-wider text-white font-serif flex items-center gap-1.5 leading-none">
              ELSEPA <span className="text-[#c5a059] font-normal italic">Simulation Suite</span>
            </h1>
            <p className="text-[11px] text-[#94a3b8] font-medium mt-1 leading-none">
              Dirac Relativistic Partial-Wave Atom Elastic Scatter & Multi-format Data Plotter
            </p>
          </div>
        </div>

        {/* Dynamic Tab Navigation HUD */}
        <nav className="flex bg-[#090a0f] border border-[#222430] p-1 rounded-xl text-xs gap-0.5">
          <button
            id="tab-simulator"
            onClick={() => setActiveTab(AppTab.SIMULATOR)}
            className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === AppTab.SIMULATOR
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <FlaskConical className="w-3.5 h-3.5 text-[#c5a059]" />
            Simulations Workstation
          </button>

          <button
            id="tab-plotter"
            onClick={() => setActiveTab(AppTab.PLOTTER)}
            className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === AppTab.PLOTTER
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <Upload className="w-3.5 h-3.5 text-[#c5a059]" />
            Spreadsheet Plotter
          </button>

          <button
            id="tab-packager"
            onClick={() => setActiveTab(AppTab.PACKAGER)}
            className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === AppTab.PACKAGER
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <Monitor className="w-3.5 h-3.5 text-[#c5a059]" />
            Desktop Executable Bundle
          </button>

          <button
            id="tab-guide"
            onClick={() => setActiveTab(AppTab.PHYSICS_GUIDE)}
            className={`px-3.5 py-2 font-bold rounded-lg transition-all cursor-pointer flex items-center gap-1.5 ${
              activeTab === AppTab.PHYSICS_GUIDE
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <BookOpen className="w-3.5 h-3.5 text-[#c5a059]" />
            Physics Guide
          </button>
        </nav>
      </header>

      {/* Main body layouts */}
      <main className="flex-1 p-6 flex flex-col gap-6 max-w-7xl mx-auto w-full self-center">
        {/* ACTIVE TAB 1: SIMULATOR WORKSPACES */}
        {activeTab === AppTab.SIMULATOR && (
          <div className="grid grid-cols-1 xl:grid-cols-12 gap-6 items-start">
            {/* Left Parameter Panel: Collapsible Accordion grids */}
            <div className="xl:col-span-3 flex flex-col gap-5 bg-[#12131a] border border-[#222430] rounded-2xl p-5 shadow-[0_4px_25px_rgba(0,0,0,0.5)] select-none">
              <div className="flex items-center justify-between border-b border-[#222430] pb-3">
                <h3 className="text-sm font-bold text-white font-serif tracking-wider flex items-center gap-2">
                  <Settings className="w-4.5 h-4.5 text-[#c5a059]" />
                  Physics Controls
                </h3>
                <button
                  id="reset-defaults-lnk"
                  onClick={handleResetDefaults}
                  className="text-[10px] text-[#c5a059] hover:text-white font-bold uppercase tracking-widest flex items-center gap-1 cursor-pointer transition-colors"
                >
                  <Undo className="w-3 h-3" />
                  Defaults
                </button>
              </div>

              {/* Accordion Group 1: Collision Incident Beam */}
              <div className="flex flex-col gap-4">
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  1. Incident Projectile
                </div>

                {/* Particle selector */}
                <div>
                  <label className="block text-xs font-semibold text-[#94a3b8] mb-1.5">
                    Projectile Particle Type
                  </label>
                  <div className="grid grid-cols-2 bg-[#090a0f] p-0.5 border border-[#222430] rounded-lg text-xs leading-none">
                    <button
                      id="opt-particle-electron"
                      onClick={() => setProjectile("electron")}
                      className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                        projectile === "electron" ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm" : "text-[#94a3b8] hover:text-white"
                      }`}
                    >
                      Electron (e⁻)
                    </button>
                    <button
                      id="opt-particle-positron"
                      onClick={() => setProjectile("positron")}
                      className={`py-1.5 rounded-md font-bold cursor-pointer transition ${
                        projectile === "positron" ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-sm" : "text-[#94a3b8] hover:text-white"
                      }`}
                    >
                      Positron (e⁺)
                    </button>
                  </div>
                </div>

                {/* Beam energy selectors */}
                <div className="flex flex-col gap-2">
                  <div className="flex items-center justify-between">
                    <label className="block text-xs font-semibold text-[#94a3b8]">
                      Kinetic Energy
                    </label>
                    <div className="flex bg-[#090a0f] border border-[#222430] p-0.5 rounded-md text-[10px] font-extrabold text-[#94a3b8]">
                      {(["eV", "keV", "MeV"] as const).map(unit => (
                        <button
                          key={unit}
                          id={`unit-${unit}`}
                          onClick={() => {
                            setEnergyUnit(unit);
                            // Set reasonable safety bounds relative to order of magnitude range
                            if (unit === "eV") setEnergy(200);
                            if (unit === "keV") setEnergy(50);
                            if (unit === "MeV") setEnergy(2);
                          }}
                          className={`px-1.5 py-0.5 rounded cursor-pointer transition ${
                            energyUnit === unit ? "bg-[#1c1d26] text-[#c5a059] font-bold" : "text-gray-500"
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
                      className="w-16 text-center text-xs border border-[#222430] p-1 rounded font-mono font-bold bg-[#090a0f] text-white focus:border-[#c5a059] outline-none"
                    />
                  </div>
                </div>
              </div>

              {/* Accordion Group 2: Advanced Relativistic Dirac Potentials */}
              <div className="flex flex-col gap-4 border-t border-[#222430] pt-4">
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  2. Dirac Potentials
                </div>

                {/* Electrostatic potential selection */}
                <div>
                  <label className="block text-xs font-semibold text-[#94a3b8] mb-1">
                    Static Potential Model
                  </label>
                  <select
                    id="select-potential-model"
                    value={potentialModel}
                    onChange={(e) => setPotentialModel(e.target.value as any)}
                    className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-bold focus:border-[#c5a059] focus:ring-1 focus:ring-[#c5a059] outline-none cursor-pointer"
                  >
                    <option value="dirac-fock">Dirac-Fock (exact)</option>
                    <option value="hartree-fock">Hartree-Fock (SCF)</option>
                    <option value="bohr-screening">Bohr Exponential Screening</option>
                    <option value="yukawa">Yukawa Potential</option>
                  </select>
                </div>

                {/* Exchange potential selection (electrons only) */}
                {projectile === "electron" && (
                  <div>
                    <label className="block text-xs font-semibold text-[#94a3b8] mb-1">
                      Exchange Interaction
                    </label>
                    <select
                      id="select-exchange-model"
                      value={exchangeModel}
                      onChange={(e) => setExchangeModel(e.target.value as any)}
                      className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-bold focus:border-[#c5a059] focus:ring-1 focus:ring-[#c5a059] outline-none cursor-pointer"
                    >
                      <option value="none">No Exchange Forces</option>
                      <option value="furness-mccarthy">Furness-McCarthy (FM)</option>
                      <option value="riley-truhlar">Riley-Truhlar (RT Local)</option>
                    </select>
                  </div>
                )}

                {/* Nuclear model options */}
                <div>
                  <label className="block text-xs font-semibold text-[#94a3b8] mb-1">
                    Nuclear Charge Density
                  </label>
                  <select
                    id="select-nuclear-model"
                    value={nuclearModel}
                    onChange={(e) => setNuclearModel(e.target.value as any)}
                    className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-bold focus:border-[#c5a059] focus:ring-1 focus:ring-[#c5a059] outline-none cursor-pointer"
                  >
                    <option value="point">Point Charge</option>
                    <option value="uniform">Uniform Sphere</option>
                    <option value="fermi">Fermi Two-parameter distribution</option>
                  </select>
                </div>
              </div>

              {/* Accordion Group 3: Real and Imaginary Auxiliary Forces */}
              <div className="flex flex-col gap-3 border-t border-[#222430] pt-4">
                <div className="text-[11px] font-bold text-[#c5a059] uppercase tracking-wider border-l-2 border-[#c5a059] pl-2 font-mono">
                  3. Dynamic Enhancements
                </div>

                <label className="flex items-center gap-2.5 cursor-pointer text-xs font-semibold text-[#e2e8f0] leading-none">
                  <input
                    id="check-correlation"
                    type="checkbox"
                    checked={correlationPolarization}
                    onChange={(e) => setCorrelationPolarization(e.target.checked)}
                    className="w-4 h-4 rounded border-[#222430] bg-[#090a0f] text-[#c5a059] focus:ring-0 focus:ring-offset-0 accent-[#c5a059]"
                  />
                  Correlation-Polarization
                </label>

                <label className="flex items-center gap-2.5 cursor-pointer text-xs font-semibold text-[#e2e8f0] leading-none">
                  <input
                    id="check-absorption"
                    type="checkbox"
                    checked={absorptionModel}
                    onChange={(e) => setAbsorptionModel(e.target.checked)}
                    className="w-4 h-4 rounded border-[#222430] bg-[#090a0f] text-[#c5a059] focus:ring-0 focus:ring-offset-0 accent-[#c5a059]"
                  />
                  Inelastic Absorption Potential
                </label>
              </div>

              {/* Scientific constants block */}
              {simulation && (
                <div className="bg-[#090a0f] border border-[#222430] rounded-xl p-3 flex flex-col gap-1.5 text-[10.5px] leading-snug text-[#94a3b8]">
                  <div className="font-bold text-white border-b border-[#222430] pb-1 mb-1 font-mono uppercase text-[9px] tracking-widest select-none flex items-center justify-between">
                    <span>COLLISION METADATA (TABULATED):</span>
                  </div>
                  <div className="flex justify-between items-center bg-[#161720] px-2 py-1.5 rounded-lg border border-[#222430] mb-1.5 select-none text-[10px]">
                    <span className="font-bold text-white uppercase font-mono">Engine Status:</span>
                    {simLoading ? (
                      <span className="text-yellow-400 font-bold font-mono animate-pulse uppercase flex items-center gap-1.5">
                        <span className="w-1.5 h-1.5 bg-yellow-400 rounded-full animate-bounce" />
                        Running...
                      </span>
                    ) : simulationEngine === "native_fortran" ? (
                      <span className="text-emerald-400 font-bold font-mono uppercase flex items-center gap-1">
                        <span className="w-1.5 h-1.5 bg-emerald-400 rounded-full animate-pulse" />
                        Native Fortran 90
                      </span>
                    ) : (
                      <span className="text-amber-500 font-bold font-mono uppercase flex items-center gap-1" title="To run native Fortran, run as Desktop App offline using any environment equipped with gfortran compiler.">
                        JS/TS Emulated
                      </span>
                    )}
                  </div>
                  <div className="flex justify-between">
                    <span>de Broglie Wavelength:</span>
                    <strong className="text-[#c5a059] font-mono">{simulation.deBroglieWavelength.toFixed(5)} Å</strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Total Elastic Cross Section:</span>
                    <strong className="text-white font-mono truncate max-w-[130px]" title={`${simulation.totalElasticCrossSection.toExponential(4)} a0²`}>
                      {simulation.totalElasticCrossSection.toExponential(4)} a₀²
                    </strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Momentum Transfer Cross Section:</span>
                    <strong className="text-white font-mono truncate max-w-[130px]" title={`${simulation.momentumTransferCrossSection.toExponential(4)} a0²`}>
                      {simulation.momentumTransferCrossSection.toExponential(4)} a₀²
                    </strong>
                  </div>
                  <div className="flex justify-between">
                    <span>Waves Summed (Lmax):</span>
                    <strong className="text-white font-mono">{simulation.phaseShifts.length - 1}</strong>
                  </div>
                </div>
              )}
            </div>

            {/* Center + Right panels: plotting, periodic table, file views */}
            <div className="xl:col-span-9 flex flex-col gap-6">
              {/* Periodic element map picker (collapsible trigger) */}
              <div className="flex flex-col gap-2">
                <div className="flex items-center justify-between bg-[#12131a] border border-[#222430] rounded-xl px-4 py-3 select-none shadow-md">
                  <span className="text-xs font-bold text-[#e2e8f0] flex items-center gap-1.5 font-mono leading-none">
                    <ListFilter className="w-4 h-4 text-[#c5a059] shrink-0" />
                    Periodic Element Mapper
                  </span>
                  <button
                    id="table-toggle-btn"
                    onClick={() => setShowPeriodicTable(!showPeriodicTable)}
                    className="text-xs text-[#c5a059] font-bold hover:text-white cursor-pointer transition-colors"
                  >
                    {showPeriodicTable ? "Collapse Table" : "Expand Table"}
                  </button>
                </div>
                {showPeriodicTable && (
                  <PeriodicTable
                    selectedZ={selectedZ}
                    onSelectElement={(atomicNumber) => setSelectedZ(atomicNumber)}
                  />
                )}
              </div>

              {/* Grand Interactive Scientific custom chart plotter */}
              <CustomChart
                dcsData={simulation?.dcsData || []}
                elementSymbol={simulation?.element.symbol || "Au"}
                energyText={`${energy} ${energyUnit}`}
                uploadedDatasets={uploadedDatasets}
                onToggleDatasetVisibility={handleToggleDatasetVisibility}
              />

              {/* Dynamic simulated file logs output panel */}
              {simulation && (
                <div className="bg-[#12131a] border border-[#222430] rounded-2xl overflow-hidden flex flex-col h-[340px] shadow-xl">
                  <div className="bg-[#090a0f] border-b border-[#222430] px-4 py-3 flex items-center justify-between select-none shrink-0">
                    <div className="flex items-center gap-3">
                      <Terminal className="w-4.5 h-4.5 text-[#c5a059] shrink-0" />
                      <div className="font-bold text-xs text-white leading-none flex bg-[#12131a] border border-[#222430] p-0.5 rounded-lg">
                        <button
                          id="btn-log-output"
                          onClick={() => setActiveFileView("output")}
                          className={`px-3 py-1 rounded-md transition cursor-pointer text-[11px] ${
                            activeFileView === "output" ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30" : "text-[#94a3b8] hover:text-white"
                          }`}
                        >
                          elsepa.out (Tabulated Output Logs)
                        </button>
                        <button
                          id="btn-log-input"
                          onClick={() => setActiveFileView("input")}
                          className={`px-3 py-1 rounded-md transition cursor-pointer text-[11px] ${
                            activeFileView === "input" ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30" : "text-[#94a3b8] hover:text-white"
                          }`}
                        >
                          elsepa.in (Input Config Options)
                        </button>
                      </div>
                    </div>

                    <button
                      id="download-log-file-btn"
                      onClick={() => 
                        activeFileView === "output" 
                           ? triggerTextFileDownload(`elsepa_${simulation.element.symbol}_${energy}${energyUnit}.out`, simulation.elsepaOutFile)
                           : triggerTextFileDownload("elsepa.in", simulation.elsepaInFile)
                      }
                      className="text-xs bg-[#1c1d26] text-white hover:bg-[#222430] border border-[#222430] rounded-lg px-3 py-1.5 font-bold cursor-pointer transition flex items-center gap-1.5"
                    >
                      <Download className="w-3.5 h-3.5 text-[#c5a059]" />
                      Download Configuration Files
                    </button>
                  </div>

                  <div className="flex-1 p-4 bg-[#040508] font-mono text-[11px] text-[#a5b4fc] overflow-auto whitespace-pre leading-normal border-t border-black">
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
                elementSymbol={simulation?.element.symbol || "Au"}
                energyText={`${energy} ${energyUnit}`}
                uploadedDatasets={uploadedDatasets}
                onToggleDatasetVisibility={handleToggleDatasetVisibility}
              />
            </div>
          </div>
        )}

        {/* ACTIVE TAB 3: DESKTOP SINGLE PORTABLE PACKAGER */}
        {activeTab === AppTab.PACKAGER && (
          <PackagerWidget />
        )}

        {/* ACTIVE TAB 4: THEORETICAL GUIDE */}
        {activeTab === AppTab.PHYSICS_GUIDE && (
          <TheoryGuide />
        )}
      </main>

      {/* Corporate science lab footer */}
      <footer className="bg-[#12131a] border-t border-[#222430] py-5 text-center text-xs text-[#94a3b8] select-none font-sans font-medium mt-auto shrink-0 shadow-[0_-4px_30px_rgba(0,0,0,0.4)]">
        <p>© 2026 eScatter ELSEPA Online Computing Suite. All scattering scattering models solved analytically under Dirac relativistic partial-wave limits.</p>
      </footer>
    </div>
  );
}
