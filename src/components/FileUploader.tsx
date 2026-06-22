/**
 * File Uploader Component with Drag & Drop
 * Parses .xlsx, .xls, .csv, or tabular .txt elastic scattering datasets using SheetJS
 */

import React, { useState, useRef } from "react";
import { UploadedDataset } from "../types";
import { Upload, FileSpreadsheet, AlertCircle, FileText, Check, Trash2, Sliders } from "lucide-react";
import * as XLSX from "xlsx";

interface FileUploaderProps {
  onAddDataset: (dataset: UploadedDataset) => void;
  uploadedDatasets: UploadedDataset[];
  onRemoveDataset: (id: string) => void;
}

// Preset nice distinct high-contrast colors for multiple dataset plotting
const DATASET_COLORS = [
  "#e11d48", // Rose
  "#16a34a", // Green
  "#9333ea", // Purple
  "#ea580c", // Orange
  "#0891b2", // Cyan
  "#db2777", // Pink
  "#ca8a04", // Yellow-Gold
];

export default function FileUploader({
  onAddDataset,
  uploadedDatasets,
  onRemoveDataset
}: FileUploaderProps) {
  const [errorMsg, setErrorMsg] = useState<string | null>(null);
  const [isDragging, setIsDragging] = useState<boolean>(false);
  const fileInputRef = useRef<HTMLInputElement | null>(null);

  // States for column mapper wizard
  const [pendingFile, setPendingFile] = useState<{
    name: string;
    headers: string[];
    records: any[];
  } | null>(null);

  const [angleCol, setAngleCol] = useState<string>("");
  const [dcsCol, setDcsCol] = useState<string>("");
  const [shermanCol, setShermanCol] = useState<string>("");

  // Process selected file
  const processFile = (file: File) => {
    setErrorMsg(null);
    const fileName = file.name;
    const extension = fileName.split(".").pop()?.toLowerCase();

    if (!["xlsx", "xls", "csv", "txt"].includes(extension || "")) {
      setErrorMsg("Unsupported file format. Please upload a valid .xlsx, .xls, .csv, or tabular .txt file.");
      return;
    }

    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const data = e.target?.result;
        if (!data) throw new Error("Could not read file data.");

        // Read workbook using SheetJS
        const workbook = XLSX.read(data, { type: "binary" });
        const firstSheetName = workbook.SheetNames[0];
        const worksheet = workbook.Sheets[firstSheetName];

        // Convert worksheet rows to raw JSON array (header row acts as keys)
        const jsonData = XLSX.utils.sheet_to_json(worksheet, { defval: "" });

        if (jsonData.length === 0) {
          throw new Error("The uploaded spreadsheet contains no data rows.");
        }

        // Get header columns
        const headers = Object.keys(jsonData[0]);

        // Attempt smart column mapping
        let detectedAngle = "";
        let detectedDcs = "";
        let detectedSherman = "";

        headers.forEach((h) => {
          const l = h.toLowerCase();
          if (l.includes("angle") || l.includes("theta") || l.includes("deg") || l === "th" || l === "x") {
            detectedAngle = h;
          }
          if (l.includes("dcs") || l.includes("cross") || l.includes("section") || l.includes("diff") || l.includes("sec") || l === "y" || l.includes("sigma")) {
            detectedDcs = h;
          }
          if (l.includes("sherman") || l.includes("polar") || l.includes("spin") || l === "s" || l.includes("sh")) {
            detectedSherman = h;
          }
        });

        // Fallbacks if not auto-detected
        if (!detectedAngle) detectedAngle = headers[0] || "";
        if (!detectedDcs) detectedDcs = headers[1] || headers[0] || "";

        // Trigger column mapping dialogue
        setPendingFile({
          name: fileName,
          headers,
          records: jsonData
        });

        setAngleCol(detectedAngle);
        setDcsCol(detectedDcs);
        setShermanCol(detectedSherman);
      } catch (err: any) {
        console.error(err);
        setErrorMsg(`Failed parsing file: ${err.message || "Spreadsheet format has structure errors."}`);
      }
    };

    reader.onerror = () => {
      setErrorMsg("Error occurred during file reading.");
    };

    reader.readAsBinaryString(file);
  };

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      processFile(e.target.files[0]);
    }
  };

  // Drag-and-drop triggers
  const handleDragOver = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = () => {
    setIsDragging(false);
  };

  const handleDrop = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    setIsDragging(false);
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      processFile(e.dataTransfer.files[0]);
    }
  };

  // Trigger input click
  const triggerFileInput = () => {
    fileInputRef.current?.click();
  };

  // Confirm Column mapping and add dataset
  const handleConfirmMapping = () => {
    if (!pendingFile || !angleCol || !dcsCol) return;

    const nextColorIndex = uploadedDatasets.length % DATASET_COLORS.length;
    const color = DATASET_COLORS[nextColorIndex];

    const newDataset: UploadedDataset = {
      id: `ds-${Date.now()}`,
      name: pendingFile.name.replace(/\.[^/.]+$/, ""), // remove extension
      fileName: pendingFile.name,
      headers: pendingFile.headers,
      records: pendingFile.records,
      angleColumn: angleCol,
      dcsColumn: dcsCol,
      shermanColumn: shermanCol || undefined,
      color,
      visible: true
    };

    onAddDataset(newDataset);
    setPendingFile(null);
    if (fileInputRef.current) fileInputRef.current.value = "";
  };

  return (
    <div id="file-uploader-section" className="bg-[#12131a] border border-[#222430] rounded-xl p-5 shadow-[0_4px_25px_rgba(0,0,0,0.5)] flex flex-col gap-4 font-sans">
      <div>
        <h3 className="text-sm font-semibold text-white flex items-center gap-1.5 leading-none font-serif tracking-wider">
          📊 Drag & Drop Scatter Data Loader
        </h3>
        <p className="text-xs text-[#94a3b8] mt-1">
          Compare simulated Dirac curves with experimental data or previous runs. Supports CSV, XLSX, and XLS formats.
        </p>
      </div>

      {errorMsg && (
        <div className="flex items-start gap-2 p-3 bg-red-950/40 text-red-300 text-xs rounded-lg border border-red-900/40">
          <AlertCircle className="w-4 h-4 shrink-0 mt-0.5" />
          <div>{errorMsg}</div>
        </div>
      )}

      {/* Column mapping wizard dialogue */}
      {pendingFile && (
        <div className="border border-[#c5a059]/30 bg-[#1c1d26]/50 rounded-xl p-4 flex flex-col gap-3">
          <div className="flex items-center gap-2 text-xs font-bold text-[#c5a059] tracking-wide">
            <Sliders className="w-4 h-4 text-[#c5a059]" />
            Column Binder: <span className="font-mono bg-[#090a0f] text-white px-1.5 py-0.5 rounded border border-[#222430]">{pendingFile.name}</span>
          </div>
          
          <p className="text-[11px] text-gray-300 leading-normal">
            We found multiple columns in your file. Please select which columns represent the Angle (X-Axis) and Cross Section (Y-Axis) for graphing:
          </p>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 my-1">
            <div>
              <label className="block text-[10px] font-bold text-[#c5a059] uppercase tracking-wider mb-1 font-serif">
                Scattering Angle (θ)
              </label>
              <select
                id="select-angle-column"
                value={angleCol}
                onChange={(e) => setAngleCol(e.target.value)}
                className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-medium focus:ring-1 focus:ring-[#c5a059] focus:outline-none"
              >
                {pendingFile.headers.map((h) => (
                  <option key={h} value={h}>{h}</option>
                ))}
              </select>
            </div>

            <div>
              <label className="block text-[10px] font-bold text-[#c5a059] uppercase tracking-wider mb-1 font-serif">
                Differential Section (DCS)
              </label>
              <select
                id="select-dcs-column"
                value={dcsCol}
                onChange={(e) => setDcsCol(e.target.value)}
                className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-medium focus:ring-1 focus:ring-[#c5a059] focus:outline-none"
              >
                {pendingFile.headers.map((h) => (
                  <option key={h} value={h}>{h}</option>
                ))}
              </select>
            </div>

            <div>
              <label className="block text-[10px] font-bold text-[#c5a059] uppercase tracking-wider mb-1 font-serif">
                Spin polarization S(θ) [Opt]
              </label>
              <select
                id="select-sherman-column"
                value={shermanCol}
                onChange={(e) => setShermanCol(e.target.value)}
                className="w-full text-xs rounded border border-[#222430] p-1.5 bg-[#090a0f] text-white font-medium focus:ring-1 focus:ring-[#c5a059] focus:outline-none"
              >
                <option value="">-- None --</option>
                {pendingFile.headers.map((h) => (
                  <option key={h} value={h}>{h}</option>
                ))}
              </select>
            </div>
          </div>

          <div className="flex justify-end gap-2 mt-1">
            <button
              id="cancel-mapping-btn"
              onClick={() => setPendingFile(null)}
              className="px-3 py-1.5 text-xs text-gray-300 bg-[#12131a] hover:bg-[#222430] border border-[#222430] rounded-lg cursor-pointer font-medium"
            >
              Cancel
            </button>
            <button
              id="confirm-mapping-btn"
              onClick={handleConfirmMapping}
              className="px-3 py-1.5 text-xs text-[#090a0f] bg-[#c5a059] hover:bg-[#dfba73] rounded-lg cursor-pointer font-bold flex items-center gap-1.5 shadow-md"
            >
              <Check className="w-3.5 h-3.5" /> Bind & Add Plot
            </button>
          </div>
        </div>
      )}

      {/* Main Drag-and-drop landing target */}
      {!pendingFile && (
        <div
          id="dropzone"
          onDragOver={handleDragOver}
          onDragLeave={handleDragLeave}
          onDrop={handleDrop}
          onClick={triggerFileInput}
          className={`border-2 border-dashed rounded-xl p-6 flex flex-col items-center justify-center gap-2 cursor-pointer transition duration-200 text-center select-none ${
            isDragging 
              ? "border-[#c5a059] bg-[#1c1d26]/40" 
              : "border-[#222430] bg-[#090a0f]/20 hover:border-[#c5a059]/40 hover:bg-[#12131a]"
          }`}
        >
          <input
            ref={fileInputRef}
            type="file"
            onChange={handleFileChange}
            accept=".xlsx,.xls,.csv,.txt"
            className="hidden"
          />
          <div className="p-3 bg-[#090a0f] rounded-full border border-[#222430] shadow-md">
            <Upload className="w-5 h-5 text-[#c5a059] animate-pulse" />
          </div>
          <div className="text-xs font-bold text-gray-200">
            Drag files here, or click to browse
          </div>
          <p className="text-[11px] text-[#94a3b8] max-w-xs leading-normal">
            Accepts experimental results from .xlsx sheets, CSVs, or tabular data outputs. Autodetect binders are applied.
          </p>
        </div>
      )}

      {/* Tabulated active loaded files manager */}
      {uploadedDatasets.length > 0 && (
        <div className="border border-[#222430] rounded-lg overflow-hidden mt-2">
          <div className="bg-[#090a0f] border-b border-[#222430] p-2.5 text-[10px] font-bold text-[#c5a059] tracking-wider uppercase select-none font-serif">
            Active Plotted Datasets ({uploadedDatasets.length})
          </div>
          <div className="divide-y divide-[#222430] bg-[#12131a]">
            {uploadedDatasets.map((ds) => (
              <div key={ds.id} className="p-3 flex items-center justify-between gap-3 text-xs border border-[#222430]">
                <div className="flex items-center gap-2.5 min-w-0">
                  {/* Distinct indicator circle */}
                  <div
                    className="w-3.5 h-3.5 rounded-full border border-black/10 shrink-0"
                    style={{ backgroundColor: ds.color }}
                  />
                  <div className="min-w-0">
                    <div className="font-bold text-white truncate leading-tight flex items-center gap-1.5">
                      <FileSpreadsheet className="w-3.5 h-3.5 text-emerald-500 shrink-0" />
                      {ds.name}
                    </div>
                    <div className="text-[10px] text-[#94a3b8] font-mono mt-0.5 truncate uppercase">
                      Cols: θ=[{ds.angleColumn}] | Y=[{ds.dcsColumn}] | rows={ds.records.length}
                    </div>
                  </div>
                </div>

                <div className="flex items-center gap-1">
                  <button
                    id={`delete-dataset-${ds.id}`}
                    onClick={() => onRemoveDataset(ds.id)}
                    className="p-1 px-2 text-[10px] text-red-400 border border-red-950/40 hover:bg-red-950/20 hover:border-red-900 rounded-lg cursor-pointer flex items-center gap-1 font-semibold"
                    title="Remove dataset"
                  >
                    <Trash2 className="w-3 h-3" />
                    Delete
                  </button>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
