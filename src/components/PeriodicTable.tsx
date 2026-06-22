/**
 * Interactive Periodic Table Component
 * Designed with a clean, high-contrast, professional style
 */

import React, { useState } from "react";
import { ElementData } from "../types";
import { PERIODIC_TABLE } from "../utils/physicsSolver";

interface PeriodicTableProps {
  selectedZ: number;
  onSelectElement: (atomicNumber: number) => void;
}

export default function PeriodicTable({ selectedZ, onSelectElement }: PeriodicTableProps) {
  const [filterCategory, setFilterCategory] = useState<string>("all");

  // Category styles mapped to beautiful dark theme metallic layouts with high readability contrast
  const categoryColors: Record<string, string> = {
    "nonmetal": "bg-emerald-950/40 text-emerald-300 border-emerald-900/50 hover:bg-emerald-900/40",
    "noble-gas": "bg-indigo-950/40 text-indigo-300 border-indigo-900/50 hover:bg-indigo-900/40",
    "alkali-metal": "bg-red-950/40 text-red-300 border-red-900/50 hover:bg-red-900/40",
    "alkaline-earth-metal": "bg-orange-950/40 text-orange-300 border-orange-900/50 hover:bg-orange-900/40",
    "metalloid": "bg-purple-950/40 text-purple-300 border-purple-900/50 hover:bg-purple-900/40",
    "halogen": "bg-teal-950/40 text-teal-300 border-teal-900/50 hover:bg-teal-900/40",
    "post-transition-metal": "bg-blue-950/40 text-blue-300 border-blue-900/50 hover:bg-blue-900/40",
    "transition-metal": "bg-amber-950/40 text-[#c5a059] border-amber-900/50 hover:bg-amber-900/40",
    "lanthanide": "bg-pink-950/40 text-pink-300 border-pink-900/50 hover:bg-pink-900/40",
    "actinide": "bg-rose-950/40 text-rose-300 border-rose-900/50 hover:bg-rose-900/40"
  };

  const categories = [
    { id: "all", label: "All Elements", color: "bg-[#090a0f] border-[#222430]" },
    { id: "nonmetal", label: "Nonmetals", color: "bg-emerald-950/40 text-emerald-300 border-emerald-900/50" },
    { id: "noble-gas", label: "Noble Gases", color: "bg-indigo-950/40 text-indigo-300 border-indigo-900/50" },
    { id: "alkali-metal", label: "Alkali Metals", color: "bg-red-950/40 text-red-300 border-red-900/50" },
    { id: "alkaline-earth-metal", label: "Alkaline Earths", color: "bg-orange-950/40 text-orange-300 border-orange-900/50" },
    { id: "transition-metal", label: "Transition Metals", color: "bg-amber-950/40 text-[#c5a059] border-amber-900/50" },
    { id: "metalloid", label: "Metalloids", color: "bg-purple-950/40 text-purple-300 border-purple-900/50" },
    { id: "lanthanide", label: "Lanthanides", color: "bg-pink-950/40 text-pink-300 border-pink-900/50" },
    { id: "actinide", label: "Actinides", color: "bg-rose-950/40 text-rose-300 border-rose-900/50" }
  ];

  // Helper to place elements in traditional standard 18-column periodic layout
  const getGridPosition = (el: ElementData) => {
    // Row 1
    if (el.number === 1) return { gridColumnStart: 1, gridRowStart: 1 };
    if (el.number === 2) return { gridColumnStart: 18, gridRowStart: 1 };

    // Row 2 and 3
    if (el.number >= 3 && el.number <= 4) return { gridColumnStart: el.number - 2, gridRowStart: 2 };
    if (el.number >= 5 && el.number <= 10) return { gridColumnStart: el.number + 8, gridRowStart: 2 };

    if (el.number >= 11 && el.number <= 12) return { gridColumnStart: el.number - 10, gridRowStart: 3 };
    if (el.number >= 13 && el.number <= 18) return { gridColumnStart: el.number, gridRowStart: 3 };

    // Row 4 and 5
    if (el.number >= 19 && el.number <= 36) return { gridColumnStart: el.number - 18, gridRowStart: 4 };
    if (el.number >= 37 && el.number <= 54) return { gridColumnStart: el.number - 36, gridRowStart: 5 };

    // Row 6 (contains Lanthanides Z=57-71 which sit in standard row 9)
    if (el.number === 55) return { gridColumnStart: 1, gridRowStart: 6 };
    if (el.number === 56) return { gridColumnStart: 2, gridRowStart: 6 };
    if (el.number >= 57 && el.number <= 71) {
      return { gridColumnStart: el.number - 54, gridRowStart: 9 }; // lanthanide row
    }
    if (el.number >= 72 && el.number <= 86) return { gridColumnStart: el.number - 68, gridRowStart: 6 };

    // Row 7 (contains Actinides Z=89-103 which sit in standard row 10)
    if (el.number === 87) return { gridColumnStart: 1, gridRowStart: 7 };
    if (el.number === 88) return { gridColumnStart: 2, gridRowStart: 7 };
    if (el.number >= 89 && el.number <= 103) {
      return { gridColumnStart: el.number - 86, gridRowStart: 10 }; // actinide row
    }
    
    return { gridColumn: "auto" };
  };

  const activeElement = PERIODIC_TABLE.find(el => el.number === selectedZ);

  return (
    <div id="periodic-table-section" className="bg-[#12131a] border border-[#222430] rounded-xl p-5 shadow-[0_4px_25px_rgba(0,0,0,0.5)]">
      <div className="flex flex-col md:flex-row md:items-center justify-between gap-4 mb-4">
        <div>
          <h3 className="text-sm font-semibold text-white flex items-center gap-1.5 font-serif tracking-wider">
            🎨 Target Element Workspace
          </h3>
          <p className="text-xs text-[#94a3b8] font-sans mt-0.5">
            Click an element in the periodic table to run the partial-wave simulation.
          </p>
        </div>

        {/* Selected element badge */}
        {activeElement && (
          <div className="flex items-center gap-3 bg-[#090a0f] border border-[#222430] rounded-lg p-2 md:p-3 pr-4 shadow-inner">
            <div className="w-10 h-10 flex flex-col justify-center items-center rounded-md border border-[#c5a059]/40 font-mono text-center shadow-xs bg-[#12131a]">
              <span className="text-[10px] leading-none text-gray-500 font-sans">{activeElement.number}</span>
              <span className="text-base leading-none font-bold text-[#c5a059]">{activeElement.symbol}</span>
            </div>
            <div className="text-left leading-normal font-sans">
              <div className="text-sm font-bold text-white">{activeElement.name}</div>
              <div className="text-[11px] text-[#94a3b8] font-mono">
                Z={activeElement.number} | M={activeElement.mass} u
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Filter Categories badges */}
      <div className="flex flex-wrap gap-1.5 mb-5 select-none no-scrollbar overflow-x-auto pb-1">
        {categories.map((cat) => (
          <button
            key={cat.id}
            id={`filter-${cat.id}`}
            onClick={() => setFilterCategory(cat.id)}
            className={`px-2.5 py-1 text-[10px] sm:text-xs rounded-full border transition-all cursor-pointer font-sans font-semibold ${
              filterCategory === cat.id
                ? "bg-[#1c1d26] text-[#c5a059] border-[#c5a059]/40 shadow-[0_0_8px_rgba(197,160,89,0.15)]"
                : "bg-[#090a0f] text-[#94a3b8] border-[#222430] hover:bg-[#1c1d26] hover:text-white"
            }`}
          >
            {cat.label}
          </button>
        ))}
      </div>

      {/* Grid container */}
      <div className="overflow-x-auto pb-2 scrollbar-thin select-none">
        <div 
          className="grid gap-1 mb-4" 
          style={{ 
            gridTemplateColumns: "repeat(18, minmax(28px, 1fr))", 
            gridTemplateRows: "repeat(10, minmax(32px, 1fr))",
            minWidth: "620px"
          }}
        >
          {PERIODIC_TABLE.map((el) => {
            const isSelected = el.number === selectedZ;
            const style = getGridPosition(el);
            const cardBg = categoryColors[el.category] || "bg-[#12131a] hover:bg-[#1c1d26]";
            const isFiltered = filterCategory !== "all" && el.category !== filterCategory;

            return (
              <button
                key={el.number}
                id={`el-${el.symbol}`}
                onClick={() => onSelectElement(el.number)}
                style={style}
                className={`flex flex-col justify-between items-center p-0.5 sm:p-1 border rounded-md transition-all cursor-pointer h-10 sm:h-12 relative ${cardBg} ${
                  isSelected 
                    ? "ring-2 ring-[#c5a059] border-transparent shadow-[0_0_12px_rgba(197,160,89,0.25)] scale-98 z-10 font-bold" 
                    : "border-[#222430]"
                } ${isFiltered ? "opacity-20 filter grayscale" : "opacity-100"}`}
              >
                <span className="text-[7px] sm:text-[8px] text-[#94a3b8] font-mono absolute top-0.5 left-1 leading-none">
                  {el.number}
                </span>
                <span className="text-xs sm:text-sm font-bold tracking-tight text-white mt-1 sm:mt-1.5 block">
                  {el.symbol}
                </span>
                <span className="text-[7px] text-gray-500 font-sans truncate w-full px-0.5 leading-none block text-center mb-0.5">
                  {el.name}
                </span>
              </button>
            );
          })}

          {/* Table Spacer for visual reference labels */}
          <div style={{ gridColumn: "3 / 13", gridRow: "1 / 4", display: "flex", justifySelf: "center", alignSelf: "center" }} className="select-none pointer-events-none text-center">
            <div className="text-[11px] font-sans text-[#c5a059]/20 font-bold tracking-widest uppercase pl-4">
              eScatter ELSEPA Engine
            </div>
          </div>

          {/* Lanthanide Spacer Row Label */}
          <div style={{ gridColumn: "1 / 3", gridRowStart: 9, display: "flex", alignItems: "center", justifyContent: "end" }} className="pr-2 select-none text-[8px] font-mono text-pink-500 uppercase font-medium">
            57-71 *
          </div>

          {/* Actinide Spacer Row Label */}
          <div style={{ gridColumn: "1 / 3", gridRowStart: 10, display: "flex", alignItems: "center", justifyContent: "end" }} className="pr-2 select-none text-[8px] font-mono text-rose-500 uppercase font-medium">
            89-103 **
          </div>
        </div>
      </div>
    </div>
  );
}
