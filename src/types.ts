/**
 * Types & Interfaces for ELSEPA Simulation Suite
 */

export enum AppTab {
  SIMULATOR = "simulator",
  PLOTTER = "plotter",
  PACKAGER = "packager",
  PHYSICS_GUIDE = "guide",
}

export type ProjectileType = "electron" | "positron";

export interface ElementData {
  number: number; // Atomic number Z
  symbol: string;
  name: string;
  mass: number;
  category: string;
  electroNegativity?: number;
  electronConfig?: string;
  group: number;
  period: number;
  block: string;
}

export type PotentialModelType = 
  | "dirac-fock" 
  | "hartree-fock" 
  | "bohr-screening" 
  | "yukawa";

export type NuclearModelType = 
  | "point" 
  | "uniform" 
  | "fermi";

export type ExchangeModelType = 
  | "none" 
  | "furness-mccarthy" 
  | "riley-truhlar";

export interface SimulationParams {
  atomicNumber: number;
  projectile: ProjectileType;
  energy: number; // in numerical value
  energyUnit: "eV" | "keV" | "MeV";
  potentialModel: PotentialModelType;
  nuclearModel: NuclearModelType;
  exchangeModel: ExchangeModelType;
  absorptionModel: boolean;
  correlationPolarization: boolean;
  gridPoints: number;
  mode?: "single" | "compound";
  compoundFormula?: string;
  compoundName?: string;
  compoundAtoms?: { atomicNumber: number; stoichiometry: number }[];
}

export interface DataPoint {
  angle: number; // Theta from 0 to 180 deg
  dcs: number; // Differential Cross Section in a0^2/sr
  dcsRutherford: number; // Standard Rutherford DCS for comparison
  Sherman: number; // Sherman spin function S(theta) from -1 to 1
}

export interface PhaseShiftPoint {
  l: number; // orbital angular momentum
  delta: number; // partial wave phase shift (spin-up)
  eta: number; // partial wave phase shift (spin-down)
}

export interface SimulationResult {
  params: SimulationParams;
  element: ElementData; // For compounds, this can represent a virtual element descriptor
  dcsData: DataPoint[];
  phaseShifts: PhaseShiftPoint[];
  totalElasticCrossSection: number; // in a0^2 or cm^2
  momentumTransferCrossSection: number; // in a0^2 or cm^2
  deBroglieWavelength: number; // in Angstroms
  elsepaInFile: string;
  elsepaOutFile: string;
  timestamp: string;
}

export interface UploadedDataset {
  id: string;
  name: string;
  fileName: string;
  headers: string[];
  records: any[];
  angleColumn: string;
  dcsColumn: string;
  shermanColumn?: string;
  color: string;
  visible: boolean;
}

export interface PresetCompound {
  name: string;
  formula: string;
  atoms: { symbol: string; atomicNumber: number; stoichiometry: number }[];
}

export interface SavedProfile {
  id: string;
  name: string;
  color: string;
  visible: boolean;
  params: SimulationParams;
  result: SimulationResult;
}
