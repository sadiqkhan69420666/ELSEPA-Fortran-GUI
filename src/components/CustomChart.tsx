/**
 * High-performance Interactive Physics Plotter
 * Custom SVG engine featuring clickable box-zoom, dual-axes, linear/log toggles,
 * and high-resolution PNG/SVG downloads.
 */

import React, { useState, useRef, useEffect, useMemo } from "react";
import { DataPoint, UploadedDataset, SavedProfile } from "../types";
import { Download, ZoomIn, ZoomOut, RotateCcw, LineChart, Shield, Sparkles } from "lucide-react";

interface CustomChartProps {
  dcsData: DataPoint[]; // simulated DCS data
  elementSymbol: string;
  energyText: string;
  uploadedDatasets: UploadedDataset[];
  onToggleDatasetVisibility: (id: string) => void;
  savedProfiles?: SavedProfile[];
  onToggleProfileVisibility?: (id: string) => void;
  onRemoveProfile?: (id: string) => void;
}

export default function CustomChart({
  dcsData,
  elementSymbol,
  energyText,
  uploadedDatasets,
  onToggleDatasetVisibility,
  savedProfiles = [],
  onToggleProfileVisibility = () => {},
  onRemoveProfile = () => {}
}: CustomChartProps) {
  // Chart Tabs: DCS (Differential Cross Section) vs Sherman (Spin Sherman Coefficient)
  const [activeMetric, setActiveMetric] = useState<"dcs" | "sherman">("dcs");
  const [yScaleType, setYScaleType] = useState<"log" | "linear">("log");

  // Plot bounds state (zoomed region)
  const [zoomRangeX, setZoomRangeX] = useState<[number, number]>([0, 180]);
  const [zoomRangeY, setZoomRangeY] = useState<[number, number] | null>(null);

  // Screened Rutherford visibility toggle
  const [showRutherford, setShowRutherford] = useState<boolean>(true);
  const [showSimulated, setShowSimulated] = useState<boolean>(true);

  // References for drag-to-zoom box selection
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState<{ x: number; y: number } | null>(null);
  const [dragEnd, setDragEnd] = useState<{ x: number; y: number } | null>(null);

  // Active hover tracking
  const [hoverPosition, setHoverPosition] = useState<{ x: number; y: number } | null>(null);
  const [hoverData, setHoverData] = useState<{
    angle: number;
    simulatedDcs?: number;
    rutherfordDcs?: number;
    simulatedSherman?: number;
    uploads: { name: string; val: number; color: string }[];
    profiles: { name: string; val: number; color: string }[];
  } | null>(null);

  // Responsive padding
  const padding = { top: 40, right: 120, bottom: 50, left: 65 };
  const graphWidth = 720;
  const graphHeight = 400;
  const viewWidth = graphWidth + padding.left + padding.right;
  const viewHeight = graphHeight + padding.top + padding.bottom;

  // Determine active datasets to draw
  const activeUploads = useMemo(() => {
    return uploadedDatasets.filter(ds => ds.visible);
  }, [uploadedDatasets]);

  // Compute total boundaries across simulated and uploaded data
  const dataBoundaries = useMemo(() => {
    // X boundaries are strictly 0 to 180 degrees
    let minX = 0;
    let maxX = 180;

    // Standard starting Y boundaries for DCS
    let minY = 1e-4;
    let maxY = 1e4;

    if (activeMetric === "dcs") {
      // Find min & max across active lines
      const values: number[] = [];
      if (showSimulated && dcsData.length > 0) {
        dcsData.forEach(d => {
          values.push(d.dcs);
          if (showRutherford) values.push(d.dcsRutherford);
        });
      }

      // Add uploaded dataset boundaries
      activeUploads.forEach(ds => {
        ds.records.forEach(r => {
          const val = parseFloat(r[ds.dcsColumn]);
          if (!isNaN(val) && val > 0) values.push(val);
        });
      });

      // Add saved overlay profiles boundaries
      savedProfiles.filter(p => p.visible).forEach(p => {
        p.result.dcsData.forEach(d => {
          const val = d.dcs;
          if (val > 0) values.push(val);
        });
      });

      if (values.length > 0) {
        minY = Math.min(...values);
        maxY = Math.max(...values);
        
        // Add padding in logarithmic or linear spaces
        if (yScaleType === "log") {
          minY = Math.pow(10, Math.floor(Math.log10(minY)) - 0.2);
          maxY = Math.pow(10, Math.ceil(Math.log10(maxY)) + 0.2);
          // Clamp lower boundary to avoid infinite logs
          if (minY < 1e-10) minY = 1e-10;
        } else {
          const span = maxY - minY;
          minY = Math.max(0, minY - span * 0.05);
          maxY = maxY + span * 0.05;
        }
      }
    } else {
      // Sherman is strictly restricted from -1 to 1 space
      minY = -1.05;
      maxY = 1.05;
    }

    return { minX, maxX, minY, maxY };
  }, [dcsData, activeMetric, showRutherford, showSimulated, yScaleType, activeUploads]);

  // Handle current dynamic zoom boundaries
  const currentXBounds = useMemo(() => zoomRangeX, [zoomRangeX]);
  const currentYBounds = useMemo(() => {
    if (zoomRangeY) return zoomRangeY;
    return [dataBoundaries.minY, dataBoundaries.maxY] as [number, number];
  }, [zoomRangeY, dataBoundaries]);

  // Reset zoom action
  const handleResetZoom = () => {
    setZoomRangeX([0, 180]);
    setZoomRangeY(null);
  };

  // Quick Preset Actions
  const handleXPreset = (min: number, max: number) => {
    setZoomRangeX([min, max]);
  };

  // Convert angular data coordinate to SVG container coordinate
  const getX = (theta: number) => {
    const minTheta = currentXBounds[0];
    const maxTheta = currentXBounds[1];
    return padding.left + ((theta - minTheta) / (maxTheta - minTheta)) * graphWidth;
  };

  // Convert DCS/Y coordinate to SVG container coordinate
  const getY = (yVal: number) => {
    const minY = currentYBounds[0];
    const maxY = currentYBounds[1];

    if (activeMetric === "dcs" && yScaleType === "log") {
      // Safeguard negative logs
      const safeVal = Math.max(minY, yVal);
      const logMin = Math.log10(minY);
      const logMax = Math.log10(maxY);
      const logVal = Math.log10(safeVal);
      return padding.top + (1 - (logVal - logMin) / (logMax - logMin)) * graphHeight;
    } else {
      return padding.top + (1 - (yVal - minY) / (maxY - minY)) * graphHeight;
    }
  };

  // Convert SVG pixel coordinate back to Angular physics coordinate
  const getThetaValue = (svgX: number) => {
    const minTheta = currentXBounds[0];
    const maxTheta = currentXBounds[1];
    const pct = (svgX - padding.left) / graphWidth;
    const value = minTheta + pct * (maxTheta - minTheta);
    return Math.max(0, Math.min(180, value));
  };

  // Convert SVG pixel coordinate back to DCS physics value
  const getYValue = (svgY: number) => {
    const minY = currentYBounds[0];
    const maxY = currentYBounds[1];
    const pct = 1 - (svgY - padding.top) / graphHeight;

    if (activeMetric === "dcs" && yScaleType === "log") {
      const logMin = Math.log10(minY);
      const logMax = Math.log10(maxY);
      return Math.pow(10, logMin + pct * (logMax - logMin));
    } else {
      return minY + pct * (maxY - minY);
    }
  };

  // Generate ticks for Gridlines
  const xTicks = useMemo(() => {
    const [min, max] = currentXBounds;
    const range = max - min;
    let step = 30;
    if (range < 15) step = 2;
    else if (range < 30) step = 5;
    else if (range < 60) step = 10;
    else if (range < 100) step = 20;

    const first = Math.ceil(min / step) * step;
    const ticks: number[] = [];
    for (let t = first; t <= max; t += step) {
      ticks.push(parseFloat(t.toFixed(2)));
    }
    return ticks;
  }, [currentXBounds]);

  const yTicks = useMemo(() => {
    const [min, max] = currentYBounds;
    const ticks: number[] = [];

    if (activeMetric === "dcs" && yScaleType === "log") {
      const firstExp = Math.ceil(Math.log10(min));
      const lastExp = Math.floor(Math.log10(max));
      for (let exp = firstExp; exp <= lastExp; exp++) {
        ticks.push(Math.pow(10, exp));
      }
      // If range is very small, generate fractional subdivisions
      if (ticks.length <= 2) {
        ticks.push(min);
        ticks.push(max);
      }
    } else {
      // Linear ticks
      const span = max - min;
      let step = span / 5;
      if (activeMetric === "sherman") step = 0.5;

      const orderOfMag = Math.pow(10, Math.floor(Math.log10(step)));
      const normalizedStep = Math.round(step / orderOfMag) * orderOfMag;
      
      const start = Math.ceil(min / normalizedStep) * normalizedStep;
      for (let t = start; t <= max; t += normalizedStep) {
        ticks.push(parseFloat(t.toFixed(4)));
      }
    }
    return ticks;
  }, [currentYBounds, yScaleType, activeMetric]);

  // Compute simulated SVG lines
  const simulatedPath = useMemo(() => {
    if (!showSimulated || dcsData.length === 0) return "";
    let path = "";
    
    dcsData.forEach((d) => {
      const sx = getX(d.angle);
      const val = activeMetric === "dcs" ? d.dcs : d.Sherman;
      const sy = getY(val);

      if (d.angle === 0) {
        path = `M ${sx} ${sy}`;
      } else {
        path += ` L ${sx} ${sy}`;
      }
    });
    return path;
  }, [dcsData, showSimulated, activeMetric, yScaleType, currentXBounds, currentYBounds]);

  const rutherfordPath = useMemo(() => {
    if (activeMetric !== "dcs" || !showRutherford || dcsData.length === 0) return "";
    let path = "";

    dcsData.forEach((d) => {
      const sx = getX(d.angle);
      const sy = getY(d.dcsRutherford);

      if (d.angle === 0) {
        path = `M ${sx} ${sy}`;
      } else {
        path += ` L ${sx} ${sy}`;
      }
    });
    return path;
  }, [dcsData, showRutherford, activeMetric, yScaleType, currentXBounds, currentYBounds]);

  // Generate SVG path coordinate line for each uploaded dataset
  const uploadedPaths = useMemo(() => {
    return activeUploads.map((ds) => {
      let path = "";
      // Sort records by angle strictly before plotting
      const sortedRecords = [...ds.records].sort((a, b) => {
        return parseFloat(a[ds.angleColumn]) - parseFloat(b[ds.angleColumn]);
      });

      let pointCount = 0;
      sortedRecords.forEach((r) => {
        const theta = parseFloat(r[ds.angleColumn]);
        const dcsRaw = parseFloat(r[ds.dcsColumn]);
        
        let val = dcsRaw;
        if (activeMetric === "sherman" && ds.shermanColumn) {
          val = parseFloat(r[ds.shermanColumn]) || 0;
        }

        if (isNaN(theta) || isNaN(val)) return;

        const sx = getX(theta);
        const sy = getY(val);

        if (pointCount === 0) {
          path = `M ${sx} ${sy}`;
        } else {
          path += ` L ${sx} ${sy}`;
        }
        pointCount++;
      });

      return {
        id: ds.id,
        name: ds.name,
        color: ds.color,
        path
      };
    });
  }, [activeUploads, activeMetric, yScaleType, currentXBounds, currentYBounds]);

  // Generate SVG path coordinate line for each overlay saved profile
  const profilePaths = useMemo(() => {
    return savedProfiles.filter(p => p.visible).map((p) => {
      let path = "";
      let pointCount = 0;
      p.result.dcsData.forEach((d) => {
        const val = activeMetric === "dcs" ? d.dcs : d.Sherman;
        const sx = getX(d.angle);
        const sy = getY(val);

        if (pointCount === 0) {
          path = `M ${sx} ${sy}`;
        } else {
          path += ` L ${sx} ${sy}`;
        }
        pointCount++;
      });
      return {
        id: p.id,
        name: p.name,
        color: p.color,
        path
      };
    });
  }, [savedProfiles, activeMetric, yScaleType, currentXBounds, currentYBounds]);

  // Handle Drag Selection Box (Box Zooming)
  const handleMouseDown = (e: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current) return;
    const rect = svgRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    // Verify pointer clicks inside the active plotting canvas box
    if (x >= padding.left && x <= padding.left + graphWidth &&
        y >= padding.top && y <= padding.top + graphHeight) {
      setIsDragging(true);
      setDragStart({ x, y });
      setDragEnd({ x, y });
    }
  };

  const handleMouseMove = (e: React.MouseEvent<SVGSVGElement>) => {
    if (!svgRef.current) return;
    const rect = svgRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    // Track active selection dragging
    if (isDragging && dragStart) {
      // Clamp values within chart viewport box limits
      const cx = Math.max(padding.left, Math.min(padding.left + graphWidth, x));
      const cy = Math.max(padding.top, Math.min(padding.top + graphHeight, y));
      setDragEnd({ x: cx, y: cy });
    }

    // Capture coords for Hover Inspection Tooltip
    if (x >= padding.left && x <= padding.left + graphWidth &&
        y >= padding.top && y <= padding.top + graphHeight) {
      
      const hoverTheta = getThetaValue(x);
      
      // Interpolate simulated entry
      const closestSimulated = dcsData.reduce((prev, curr) => {
        return Math.abs(curr.angle - hoverTheta) < Math.abs(prev.angle - hoverTheta) ? curr : prev;
      }, dcsData[0] || { angle: 0, dcs: 1, dcsRutherford: 1, Sherman: 0 });

      // Find closest values in uploaded datasets too!
      const hoverUploads = activeUploads.map(ds => {
        let closestRecordVal = 0;
        let diff = Infinity;
        ds.records.forEach(r => {
          const t = parseFloat(r[ds.angleColumn]);
          const d_val = parseFloat(r[activeMetric === "dcs" ? ds.dcsColumn : (ds.shermanColumn || "")]);
          if (!isNaN(t) && !isNaN(d_val)) {
            const currentDiff = Math.abs(t - hoverTheta);
            if (currentDiff < diff) {
              diff = currentDiff;
              closestRecordVal = d_val;
            }
          }
        });
        return { name: ds.name, val: closestRecordVal, color: ds.color };
      });

      // Find closest values in overlaid saved profiles too!
      const hoverProfiles = savedProfiles.filter(p => p.visible).map(p => {
        let closestPt = p.result.dcsData.reduce((prev, curr) => {
          return Math.abs(curr.angle - hoverTheta) < Math.abs(prev.angle - hoverTheta) ? curr : prev;
        }, p.result.dcsData[0] || { angle: 0, dcs: 1, dcsRutherford: 1, Sherman: 0 });
        const val = activeMetric === "dcs" ? closestPt.dcs : closestPt.Sherman;
        return { name: p.name, val, color: p.color };
      });

      setHoverPosition({ x, y });
      setHoverData({
        angle: parseFloat(hoverTheta.toFixed(2)),
        simulatedDcs: closestSimulated.dcs,
        rutherfordDcs: closestSimulated.dcsRutherford,
        simulatedSherman: closestSimulated.Sherman,
        uploads: hoverUploads,
        profiles: hoverProfiles
      });
    } else {
      setHoverPosition(null);
      setHoverData(null);
    }
  };

  const handleMouseUp = () => {
    if (isDragging && dragStart && dragEnd) {
      setIsDragging(false);

      // Verify dragged viewport is significant size (at least 5px drag grid delta)
      const dx = Math.abs(dragEnd.x - dragStart.x);
      const dy = Math.abs(dragEnd.y - dragStart.y);

      if (dx > 5 && dy > 5) {
        // Compute new numerical ranges for axes
        const theta1 = getThetaValue(Math.min(dragStart.x, dragEnd.x));
        const theta2 = getThetaValue(Math.max(dragStart.x, dragEnd.x));

        const y1 = getYValue(Math.max(dragStart.y, dragEnd.y));
        const y2 = getYValue(Math.min(dragStart.y, dragEnd.y));

        setZoomRangeX([theta1, theta2]);
        setZoomRangeY([y1, y2]);
      }
      setDragStart(null);
      setDragEnd(null);
    }
  };

  // Trigger downloads (Vector SVG or high DPI PNG from offscreen canvas)
  const downloadChartAsVectorSvg = () => {
    if (!svgRef.current) return;
    
    // Create static standalone copy with explicit styles
    const svgContent = svgRef.current.outerHTML;
    const blob = new Blob([
      `<?xml version="1.0" standalone="no"?>\n<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n${svgContent}`
    ], { type: "image/svg+xml;charset=utf-8" });

    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `ELSEPA_Plot_${elementSymbol}_${energyText}_${activeMetric}.svg`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const downloadChartAsHighResPng = () => {
    if (!svgRef.current) return;

    // Wait, let's render the detailed plot onto an HTML canvas using CanvasRenderingContext2D.
    // Highly accurate canvas rendering generates extremely polished, crisp corporate science report images!
    const canvas = document.createElement("canvas");
    // Generate high resolution 3x scaling for sharp text and subpixel curve drawing
    const scale = 3;
    canvas.width = viewWidth * scale;
    canvas.height = viewHeight * scale;
    
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    ctx.scale(scale, scale);
    
    // Draw white canvas background
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, viewWidth, viewHeight);

    // Draw main plot area background
    ctx.fillStyle = "#fcfdfd";
    ctx.fillRect(padding.left, padding.top, graphWidth, graphHeight);

    // Draw Axis grid lines
    ctx.strokeStyle = "#f1f3f5";
    ctx.lineWidth = 0.8;
    ctx.font = "500 10px Inter, sans-serif";
    ctx.fillStyle = "#4b5563";

    // Draw grid columns
    xTicks.forEach(t => {
      const cx = getX(t);
      ctx.beginPath();
      ctx.moveTo(cx, padding.top);
      ctx.lineTo(cx, padding.top + graphHeight);
      ctx.stroke();

      // Tick Labels
      ctx.textAlign = "center";
      ctx.fillText(`${t}°`, cx, padding.top + graphHeight + 16);
    });

    // Draw grid rows
    yTicks.forEach(t => {
      const cy = getY(t);
      ctx.beginPath();
      ctx.moveTo(padding.left, cy);
      ctx.lineTo(padding.left + graphWidth, cy);
      ctx.stroke();

      // Tick Labels
      ctx.textAlign = "right";
      const displayVal = activeMetric === "dcs" && yScaleType === "log" 
        ? t.toExponential(0) 
        : t.toString();
      ctx.fillText(displayVal, padding.left - 8, cy + 3.5);
    });

    // Draw outer frame boundaries
    ctx.strokeStyle = "#e5e7eb";
    ctx.lineWidth = 1.2;
    ctx.strokeRect(padding.left, padding.top, graphWidth, graphHeight);

    // Labels & titles
    ctx.font = "bold 13px Inter, sans-serif";
    ctx.fillStyle = "#111827";
    ctx.textAlign = "center";
    ctx.fillText(
      `ELSEPA Scattering Analysis: ${elementSymbol} at ${energyText} (${activeMetric.toUpperCase()} Plot)`,
      padding.left + graphWidth / 2,
      padding.top - 18
    );

    // Axis titles
    ctx.font = "500 11px Inter, sans-serif";
    ctx.fillText("Scattering Angle (θ)", padding.left + graphWidth / 2, padding.top + graphHeight + 38);

    ctx.save();
    ctx.translate(20, padding.top + graphHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(
      activeMetric === "dcs" 
        ? "Differential Cross Section (a0^2/sr)" 
        : "Sherman Spin polarization function S(θ)",
      0, 0
    );
    ctx.restore();

    // FUNCTION helper to plot SVG path to canvas
    const drawSvgPathOnCanvas = (pathString: string, strokeColor: string, strokeWidth: number, dashed: boolean = false) => {
      if (!pathString) return;
      ctx.save();
      ctx.strokeStyle = strokeColor;
      ctx.lineWidth = strokeWidth;
      if (dashed) ctx.setLineDash([5, 4]);

      ctx.beginPath();
      const tokens = pathString.split(/\s+/);
      let currentX = 0;
      let currentY = 0;

      for (let i = 0; i < tokens.length; i++) {
        const tok = tokens[i];
        if (tok === "M" || tok === "L") {
          const targetX = parseFloat(tokens[i + 1]);
          const targetY = parseFloat(tokens[i + 2]);
          if (!isNaN(targetX) && !isNaN(targetY)) {
            // Clip coordinates inside chart viewport bounding box physically
            const boundedX = Math.max(padding.left, Math.min(padding.left + graphWidth, targetX));
            const boundedY = Math.max(padding.top, Math.min(padding.top + graphHeight, targetY));
            if (tok === "M") {
              ctx.moveTo(boundedX, boundedY);
            } else {
              ctx.lineTo(boundedX, boundedY);
            }
          }
          i += 2;
        }
      }
      ctx.stroke();
      ctx.restore();
    };

    // Plot Screened Rutherford if flagged
    if (activeMetric === "dcs" && showRutherford) {
      drawSvgPathOnCanvas(rutherfordPath, "#9ca3af", 1.5, true);
    }

    // Plot active simulated results
    if (showSimulated) {
      drawSvgPathOnCanvas(simulatedPath, "#2563eb", 2.2, false);
    }

    // Plot all active uploads
    uploadedPaths.forEach(up => {
      drawSvgPathOnCanvas(up.path, up.color, 2.0, false);
    });

    // Draw simple crisp publication legends
    let legendX = padding.left + graphWidth + 12;
    let legendY = padding.top + 20;
    ctx.font = "600 10px Inter, sans-serif";
    ctx.fillStyle = "#374151";
    ctx.textAlign = "left";
    ctx.fillText("DATASETS:", legendX, legendY);

    legendY += 16;
    if (showSimulated) {
      ctx.fillStyle = "#2563eb";
      ctx.fillRect(legendX, legendY - 6, 12, 1.8);
      ctx.fillStyle = "#4b5563";
      ctx.font = "500 9px Inter, sans-serif";
      ctx.fillText(`Simulated (${elementSymbol})`, legendX + 18, legendY);
      legendY += 16;
    }

    if (activeMetric === "dcs" && showRutherford) {
      ctx.strokeStyle = "#9ca3af";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(legendX, legendY - 5);
      ctx.lineTo(legendX + 12, legendY - 5);
      ctx.stroke();
      ctx.fillStyle = "#4b5563";
      ctx.fillText("Rutherford (Screened)", legendX + 18, legendY);
      legendY += 16;
    }

    uploadedPaths.forEach(up => {
      ctx.fillStyle = up.color;
      ctx.fillRect(legendX, legendY - 6, 12, 2.0);
      ctx.fillStyle = "#4b5563";
      ctx.fillText(up.name, legendX + 18, legendY);
      legendY += 16;
    });

    // Generate link and download image file
    const url = canvas.toDataURL("image/png");
    const link = document.createElement("a");
    link.href = url;
    link.download = `ELSEPA_Publication_Plot_${elementSymbol}_${energyText}.png`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <div id="physics-analytics-plot-card" className="bg-[#12131a] border border-[#222430] rounded-xl p-5 shadow-[0_4px_25px_rgba(0,0,0,0.5)] flex flex-col gap-4">
      {/* Tab controls */}
      <div className="flex flex-col sm:flex-row sm:items-center justify-between border-b border-[#222430] pb-3 gap-3">
        <div className="flex bg-[#090a0f] border border-[#222430] p-0.5 rounded-lg select-none">
          <button
            id="tab-metric-dcs"
            onClick={() => { setActiveMetric("dcs"); handleResetZoom(); }}
            className={`px-3 py-1.5 text-xs font-semibold rounded-md transition-all cursor-pointer flex items-center gap-1.5 ${
              activeMetric === "dcs"
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-[0_0_8px_rgba(197,160,89,0.15)]"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <LineChart className="w-3.5 h-3.5 text-[#c5a059]" />
            Differential Cross Section (DCS)
          </button>
          <button
            id="tab-metric-sherman"
            onClick={() => { setActiveMetric("sherman"); handleResetZoom(); }}
            className={`px-3 py-1.5 text-xs font-semibold rounded-md transition-all cursor-pointer flex items-center gap-1.5 ${
              activeMetric === "sherman"
                ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/30 shadow-[0_0_8px_rgba(197,160,89,0.15)]"
                : "text-[#94a3b8] hover:text-white"
            }`}
          >
            <Shield className="w-3.5 h-3.5 text-amber-500" />
            Sherman S(θ) Polarization
          </button>
        </div>

        {/* Dynamic plotting scale selection info */}
        <div className="flex items-center gap-2">
          {activeMetric === "dcs" && (
            <div className="flex bg-[#090a0f] border border-[#222430] p-0.5 rounded-lg text-xs font-semibold select-none">
              <button
                id="toggle-scale-log"
                onClick={() => setYScaleType("log")}
                className={`px-2 py-1 rounded-md transition-all cursor-pointer ${
                  yScaleType === "log" ? "bg-[#1b1c25] text-[#c5a059] font-bold border border-[#c5a059]/20" : "text-[#94a3b8] hover:text-white"
                }`}
              >
                Logscale
              </button>
              <button
                id="toggle-scale-linear"
                onClick={() => setYScaleType("linear")}
                className={`px-2 py-1 rounded-md transition-all cursor-pointer ${
                  yScaleType === "linear" ? "bg-[#1b1c25] text-[#c5a059] font-bold border border-[#c5a059]/20" : "text-[#94a3b8] hover:text-white"
                }`}
              >
                Linear
              </button>
            </div>
          )}

          {/* Quick interactive zoom actions */}
          <button
            id="zoom-reset-btn"
            onClick={handleResetZoom}
            title="Reset range"
            className="flex items-center gap-1.5 px-3 py-1.5 text-xs text-white border border-[#222430] rounded-lg bg-[#1c1d26] hover:bg-[#222430] cursor-pointer"
          >
            <RotateCcw className="w-3.5 h-3.5 text-[#c5a059]" />
            Reset Zoom
          </button>
        </div>
      </div>

      {/* Quick Angular range zooms */}
      <div className="flex flex-wrap items-center gap-1.5 text-xs select-none">
        <span className="text-[#94a3b8] font-sans">Quick Angle Ranges:</span>
        <button
          onClick={() => handleXPreset(0, 180)}
          className="px-2.5 py-1 rounded border border-[#222430] bg-[#1c1d26] hover:bg-[#222430] text-gray-300 cursor-pointer text-[11px]"
        >
          Full (0° - 180°)
        </button>
        <button
          onClick={() => handleXPreset(0, 15)}
          className="px-2.5 py-1 rounded border border-[#222430] bg-[#1c1d26] hover:bg-[#222430] text-gray-300 cursor-pointer text-[11px] flex items-center gap-1"
        >
          <ZoomIn className="w-2.5 h-2.5 text-[#c5a059]" /> Forward Res (0° - 15°)
        </button>
        <button
          onClick={() => handleXPreset(15, 90)}
          className="px-2.5 py-1 rounded border border-[#222430] bg-[#1c1d26] hover:bg-[#222430] text-gray-300 cursor-pointer text-[11px]"
        >
          Mid (15° - 90°)
        </button>
        <button
          onClick={() => handleXPreset(90, 180)}
          className="px-2.5 py-1 rounded border border-[#222430] bg-[#1c1d26] hover:bg-[#222430] text-gray-300 cursor-pointer text-[11px]"
        >
          Backscattering (90° - 180°)
        </button>
      </div>

      {/* The main SVG Plotting Window */}
      <div className="relative border border-[#222430] rounded-lg bg-[#090a0f] overflow-hidden flex justify-center py-2.5 select-none shadow-inner">
        <svg
          ref={svgRef}
          viewBox={`0 0 ${viewWidth} ${viewHeight}`}
          className="w-full max-w-4xl h-auto"
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          style={{ cursor: isDragging ? "crosshair" : "default" }}
        >
          {/* Main graph grid background */}
          <rect
            x={padding.left}
            y={padding.top}
            width={graphWidth}
            height={graphHeight}
            fill="#0b0c10"
            stroke="#222430"
            strokeWidth="0.8"
          />

          {/* Grid lines & ticks (X-Axis) */}
          {xTicks.map((tick) => {
            const x = getX(tick);
            return (
              <g key={`x-${tick}`}>
                <line
                  x1={x}
                  y1={padding.top}
                  x2={x}
                  y2={padding.top + graphHeight}
                  stroke="#161821"
                  strokeWidth="0.8"
                />
                <text
                  x={x}
                  y={padding.top + graphHeight + 16}
                  textAnchor="middle"
                  className="font-mono text-[10px]"
                  fill="#94a3b8"
                >
                  {tick}°
                </text>
              </g>
            );
          })}

          {/* Grid lines & ticks (Y-Axis) */}
          {yTicks.map((tick) => {
            const y = getY(tick);
            // Format labels cleanly
            const labelText = activeMetric === "dcs" && yScaleType === "log"
              ? tick.toExponential(0)
              : tick.toString();

            return (
              <g key={`y-${tick}`}>
                <line
                  x1={padding.left}
                  y1={y}
                  x2={padding.left + graphWidth}
                  y2={y}
                  stroke="#161821"
                  strokeWidth="0.8"
                />
                <text
                  x={padding.left - 8}
                  y={y + 3.5}
                  textAnchor="end"
                  className="font-mono text-[10px]"
                  fill="#94a3b8"
                >
                  {labelText}
                </text>
              </g>
            );
          })}

          {/* Screened Rutherford scattering path (Elegant Copper-gold Dashed) */}
          {activeMetric === "dcs" && showRutherford && rutherfordPath && (
            <path
              d={rutherfordPath}
              fill="none"
              stroke="#a08455"
              strokeWidth="1.5"
              strokeDasharray="5,4"
              className="transition-all"
            />
          )}

          {/* Tabulated simulated physical core path (Luxury Gold Solid) */}
          {showSimulated && simulatedPath && (
            <path
              d={simulatedPath}
              fill="none"
              stroke="#c5a059"
              strokeWidth="2.5"
              strokeLinecap="round"
              strokeLinejoin="round"
              className="transition-all"
            />
          )}

          {/* Custom Uploaded datasets paths */}
          {uploadedPaths.map((up) => (
            up.path && (
              <path
                key={up.id}
                d={up.path}
                fill="none"
                stroke={up.color}
                strokeWidth="2.0"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            )
          ))}

          {/* Overlaid Saved Profiles paths */}
          {profilePaths.map((pp) => (
            pp.path && (
              <path
                key={pp.id}
                d={pp.path}
                fill="none"
                stroke={pp.color}
                strokeWidth="2.0"
                strokeLinecap="round"
                strokeLinejoin="round"
                className="transition-all"
              />
            )
          ))}

          {/* Visual Overlay selection dragging box (Box Zooming bounds UI) */}
          {isDragging && dragStart && dragEnd && (
            <rect
              x={Math.min(dragStart.x, dragEnd.x)}
              y={Math.min(dragStart.y, dragEnd.y)}
              width={Math.abs(dragEnd.x - dragStart.x)}
              height={Math.abs(dragEnd.y - dragStart.y)}
              fill="rgba(197, 160, 89, 0.15)"
              stroke="#c5a059"
              strokeWidth="1"
              strokeDasharray="3,2"
            />
          )}

          {/* Vertically sliding inspection line on hovering */}
          {hoverPosition && (
            <line
              x1={hoverPosition.x}
              y1={padding.top}
              x2={hoverPosition.x}
              y2={padding.top + graphHeight}
              stroke="#c5a059"
              strokeWidth="1.2"
              strokeDasharray="4,2"
              opacity="0.6"
            />
          )}

          {/* Axis Labels */}
          <text
            x={padding.left + graphWidth / 2}
            y={padding.top + graphHeight + 38}
            textAnchor="middle"
            className="font-serif text-[11px] font-bold tracking-wider"
            fill="#e2e8f0"
          >
            Collision Scattering Angle θ (degrees)
          </text>

          <g transform={`translate(18, ${padding.top + graphHeight / 2}) rotate(-90)`}>
            <text
              textAnchor="middle"
              className="font-serif text-[11px] font-bold tracking-wider"
              fill="#e2e8f0"
            >
              {activeMetric === "dcs"
                ? `Differential Cross Section (a0²/sr) ${yScaleType === "log" ? "[Log Scale]" : "[Linear]"}`
                : "Sherman Spin polarization function S(θ) []"}
            </text>
          </g>

          {/* Embedded Legend layout */}
          <g transform={`translate(${padding.left + graphWidth + 12}, ${padding.top + 16})`}>
            <text x="0" y="0" className="font-serif text-[10px] font-bold tracking-wider" fill="#c5a059">LEGEND</text>

            {showSimulated && (
              <g transform="translate(0, 16)">
                <line x1="0" y1="-3" x2="14" y2="-3" stroke="#c5a059" strokeWidth="2.5" />
                <text x="20" y="1" className="font-sans text-[10px] font-semibold" fill="#e2e8f0">
                  Simulated ({elementSymbol})
                </text>
              </g>
            )}

            {activeMetric === "dcs" && showRutherford && (
              <g transform="translate(0, 32)">
                <line x1="0" y1="-3" x2="14" y2="-3" stroke="#a08455" strokeWidth="1.5" strokeDasharray="3,2" />
                <text x="20" y="1" className="font-sans text-[10px] font-semibold" fill="#94a3b8">
                  Rutherford S.
                </text>
              </g>
            )}

            {uploadedPaths.map((up, idx) => (
              <g key={up.id} transform={`translate(0, ${48 + idx * 16})`}>
                <line x1="0" y1="-3" x2="14" y2="-3" stroke={up.color} strokeWidth="2" />
                <text x="20" y="1" className="font-sans text-[10px] font-semibold truncate w-24" fill="#e2e8f0">
                  {up.name.length > 15 ? `${up.name.substring(0, 12)}...` : up.name}
                </text>
              </g>
            ))}

            {profilePaths.map((pp, idx) => (
              <g key={pp.id} transform={`translate(0, ${48 + (uploadedPaths.length + idx) * 16})`}>
                <line x1="0" y1="-3" x2="14" y2="-3" stroke={pp.color} strokeWidth="2" />
                <text x="20" y="1" className="font-sans text-[10px] font-semibold truncate w-24" fill="#e2e8f0">
                  {pp.name.length > 15 ? `${pp.name.substring(0, 12)}...` : pp.name}
                </text>
              </g>
            ))}
          </g>
        </svg>

        {/* Floating Tooltip HTML HUD displaying coordinate details */}
        {hoverPosition && hoverData && (
          <div
            className="absolute z-20 bg-[#090a0f]/95 backdrop-blur-xs text-white text-[11px] p-2.5 rounded-lg shadow-[0_4px_30px_rgba(0,0,0,0.8)] border border-[#c5a059]/30 leading-snug font-mono flex flex-col gap-1 pointer-events-none"
            style={{
              left: `${(hoverPosition.x / viewWidth) * 100 + 1}%`,
              top: `${(hoverPosition.y / viewHeight) * 100 - 15}%`
            }}
          >
            <div className="font-serif font-bold text-[#c5a059] border-b border-[#222430] pb-1 flex items-center gap-1.5 tracking-wider">
              <Sparkles className="w-3 h-3 text-[#c5a059]" /> Angle θ: {hoverData.angle}°
            </div>
            {activeMetric === "dcs" ? (
              <>
                {showSimulated && (
                  <div className="flex items-center justify-between gap-4">
                    <span className="text-[#dfba73]">Sim DCS:</span>
                    <span>{hoverData.simulatedDcs?.toExponential(4)}</span>
                  </div>
                )}
                {showRutherford && (
                  <div className="flex items-center justify-between gap-4 text-gray-400">
                    <span className="text-[#a08455]">Rutherford:</span>
                    <span>{hoverData.rutherfordDcs?.toExponential(4)}</span>
                  </div>
                )}
              </>
            ) : (
              showSimulated && (
                <div className="flex items-center justify-between gap-4">
                  <span className="text-amber-500">S(θ) Spin:</span>
                  <span>{hoverData.simulatedSherman?.toFixed(5)}</span>
                </div>
              )
            )}

            {/* Custom Uploads tooltip rendering */}
            {hoverData.uploads.map((up, idx) => (
              <div key={`${up.name}-${idx}`} className="flex items-center justify-between gap-4">
                <span style={{ color: up.color }}>{up.name}:</span>
                <span>{up.val.toExponential(4)}</span>
              </div>
            ))}

            {/* Custom Profiles tooltip rendering */}
            {hoverData.profiles?.map((p, idx) => (
              <div key={`${p.name}-${idx}`} className="flex items-center justify-between gap-4">
                <span style={{ color: p.color }}>{p.name}:</span>
                <span>{p.val.toExponential(4)}</span>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Legending Filters Control Grid */}
      <div className="grid grid-cols-2 lg:grid-cols-4 border border-[#222430] rounded-lg p-3 gap-3 bg-[#090a0f] select-none">
        <label className="flex items-center gap-2 cursor-pointer text-xs font-semibold text-gray-300">
          <input
            id="checkbox-dataset-simulated"
            type="checkbox"
            checked={showSimulated}
            onChange={(e) => setShowSimulated(e.target.checked)}
            className="w-4 h-4 accent-[#c5a059] border-[#222430] rounded focus:ring-[#c5a059]"
          />
          Show Simulated DCS
        </label>

        {activeMetric === "dcs" && (
          <label className="flex items-center gap-2 cursor-pointer text-xs font-semibold text-gray-300">
            <input
              id="checkbox-dataset-rutherford"
              type="checkbox"
              checked={showRutherford}
              onChange={(e) => setShowRutherford(e.target.checked)}
              className="w-4 h-4 accent-[#c5a059] border-[#222430] rounded focus:ring-[#c5a059]"
            />
            Show Screened Rutherford
          </label>
        )}

        {/* Uploaded dataset visibilities */}
        {uploadedDatasets.map((ds) => (
          <label
            key={ds.id}
            id={`label-visibility-${ds.id}`}
            className="flex items-center gap-2 cursor-pointer text-xs font-semibold text-gray-300"
          >
            <input
              id={`checkbox-visibility-${ds.id}`}
              type="checkbox"
              checked={ds.visible}
              onChange={() => onToggleDatasetVisibility(ds.id)}
              className="w-4 h-4 border-[#222430] rounded focus:ring-[#c5a059]"
              style={{ accentColor: ds.color }}
            />
            <span className="truncate w-full block" style={{ color: ds.color }}>
              Show {ds.name}
            </span>
          </label>
        ))}

        {/* Saved overlaid overlay profile visibilities */}
        {savedProfiles.map((p) => (
          <label
            key={p.id}
            id={`label-overlay-visibility-${p.id}`}
            className="flex items-center gap-2 cursor-pointer text-xs font-semibold text-gray-300"
          >
            <input
              id={`checkbox-overlay-visibility-${p.id}`}
              type="checkbox"
              checked={p.visible}
              onChange={() => onToggleProfileVisibility(p.id)}
              className="w-4 h-4 border-[#222430] rounded focus:ring-[#c5a059]"
              style={{ accentColor: p.color }}
            />
            <span className="truncate w-full block" style={{ color: p.color }}>
              Show {p.name}
            </span>
          </label>
        ))}
      </div>

      {/* High-Resolution Scientific exports */}
      <div className="flex flex-col sm:flex-row items-center justify-between gap-3 font-sans mt-1">
        <div className="text-xs text-[#94a3b8] font-medium">
          💡 Click and drag a rectangular box inside the dark chart workspace to zoom in!
        </div>

        <div className="flex flex-wrap items-center gap-2 w-full sm:w-auto">
          <button
            id="download-png-btn"
            onClick={downloadChartAsHighResPng}
            className="flex items-center justify-center gap-2 px-3 py-1.5 text-xs font-bold text-[#090a0f] bg-[#c5a059] border border-[#c5a059] rounded-lg hover:bg-[#dfba73] hover:border-[#dfba73] transition duration-200 shadow-[0_4px_12px_rgba(197,160,89,0.2)] cursor-pointer w-full sm:w-auto"
          >
            <Download className="w-3.5 h-3.5" />
            Download Publication PNG (3x High-Res)
          </button>
          
          <button
            id="download-svg-btn"
            onClick={downloadChartAsVectorSvg}
            className="flex items-center justify-center gap-2 px-3 py-1.5 text-xs font-medium text-gray-300 border border-[#222430] rounded-lg bg-[#1c1d26] hover:bg-[#222430] transition duration-200 shadow-md cursor-pointer w-full sm:w-auto"
          >
            <Download className="w-3.5 h-3.5 text-[#c5a059]" />
            Download Vector SVG
          </button>
        </div>
      </div>
    </div>
  );
}
