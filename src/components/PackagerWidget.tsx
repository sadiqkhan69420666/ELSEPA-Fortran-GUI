/**
 * 1-Click Desktop Packager Component
 * Provides detailed, clear instructions and downloads a completely automated
 * Electron wrapper package, enabling users to generate a single portable
 * desktop executable (.exe or .app) without manual coding or terminal hassle.
 */

import React, { useState, useEffect } from "react";
import { Download, FileCode, CheckCircle, HelpCircle, Monitor, ShieldCheck, Terminal, HardDrive } from "lucide-react";

interface DesktopTemplate {
  appName: string;
  appUrl: string;
  files: {
    "package.json": string;
    "main.js": string;
    "README.md": string;
  };
}

export default function PackagerWidget() {
  const [template, setTemplate] = useState<DesktopTemplate | null>(null);
  const [activeFile, setActiveFile] = useState<"package.json" | "main.js" | "README.md">("README.md");
  const [downloading, setDownloading] = useState(false);
  const [success, setSuccess] = useState(false);

  useEffect(() => {
    // Retrieve custom assembled desktop wrapper template from backend
    fetch("/api/desktop/template")
      .then((res) => res.json())
      .then((data) => {
        setTemplate(data);
      })
      .catch((err) => {
        console.error("Failed fetching desktop template:", err);
      });
  }, []);

  // Downloads a specific file directly to the client's file system
  const triggerFileDownload = (fileName: string, content: string) => {
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

  const handleDownloadAll = () => {
    if (!template) return;
    setDownloading(true);

    // Download the files consecutively so they can assemble their folder
    setTimeout(() => {
      triggerFileDownload("package.json", template.files["package.json"]);
      triggerFileDownload("main.js", template.files["main.js"]);
      triggerFileDownload("README.md", template.files["README.md"]);
      
      setDownloading(false);
      setSuccess(true);
      setTimeout(() => setSuccess(false), 5000);
    }, 800);
  };

  return (
    <div id="desktop-packager-section" className="bg-[#12131a] border border-[#222430] rounded-xl p-6 shadow-[0_4px_25px_rgba(0,0,0,0.5)] flex flex-col gap-5 font-sans">
      <div className="flex flex-col md:flex-row md:items-center justify-between border-b border-[#222430] pb-4 gap-4 select-none">
        <div>
          <h2 className="text-base font-bold text-white flex items-center gap-2 font-serif tracking-wider">
            <Monitor className="w-5 h-5 text-[#c5a059]" />
            1-Click Desktop Executable Packager
          </h2>
          <p className="text-xs text-[#94a3b8] mt-1">
            Build your private, portable ELSEPA Desktop Application (.exe / .app) targeting Windows & macOS.
          </p>
        </div>

        <button
          id="packager-download-package-btn"
          onClick={handleDownloadAll}
          disabled={!template || downloading}
          className="flex items-center justify-center gap-2 px-4.5 py-2 text-xs font-bold text-[#090a0f] bg-[#c5a059] hover:bg-[#dfba73] disabled:opacity-50 border border-[#c5a059] rounded-xl cursor-pointer transition shadow-[0_4px_12px_rgba(197,160,89,0.2)] duration-200"
        >
          {success ? (
            <>
              <CheckCircle className="w-4 h-4 text-emerald-950" />
              Downloaded Configuration!
            </>
          ) : (
            <>
              <Download className="w-4 h-4" />
              {downloading ? "Formatting..." : "Download Desktop Wrapper Files"}
            </>
          )}
        </button>
      </div>

      {/* Intro features grid */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="border border-[#222430] bg-[#090a0f]/45 p-4 rounded-xl flex flex-col gap-2">
          <Terminal className="w-5 h-5 text-purple-400 shrink-0" />
          <h4 className="text-xs font-bold text-white font-serif tracking-wide">Zero Build Complexity</h4>
          <p className="text-[11px] text-[#94a3b8] leading-normal">
            No Fortran compilation or platform binary paths required. Electron manages the package bundling automatically.
          </p>
        </div>

        <div className="border border-[#222430] bg-[#090a0f]/45 p-4 rounded-xl flex flex-col gap-2">
          <ShieldCheck className="w-5 h-5 text-emerald-400 shrink-0" />
          <h4 className="text-xs font-bold text-white font-serif tracking-wide">Portable & Independent</h4>
          <p className="text-[11px] text-[#94a3b8] leading-normal">
            Compiles into a single offline-centric desktop shell (.exe for Windows, .app for macOS) that you can distribute anywhere.
          </p>
        </div>

        <div className="border border-[#222430] bg-[#090a0f]/45 p-4 rounded-xl flex flex-col gap-2">
          <HardDrive className="w-5 h-5 text-amber-500 shrink-0" />
          <h4 className="text-xs font-bold text-white font-serif tracking-wide">Full Cloud Run Syncing</h4>
          <p className="text-[11px] text-[#94a3b8] leading-normal">
            All physics algorithms, solvers, Excel plots, and Gemini engines update instantly inside your desktop container shell.
          </p>
        </div>
      </div>

      {/* Guide details & live files preview */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-5 mt-2">
        <div className="lg:col-span-5 flex flex-col gap-3">
          <h3 className="text-xs font-bold text-[#c5a059] uppercase tracking-wider select-none font-serif">
            🚀 3-Step Binaries Compilation
          </h3>

          <ol className="list-decimal list-inside text-xs text-gray-300 flex flex-col gap-3 leading-relaxed">
            <li>
              <strong>Download Configuration</strong>: Click the blue button above to save the preconfigured wrapper folder contents (<code className="font-mono text-[#dfba73] font-bold bg-[#090a0f] border border-[#222430] px-1 rounded">package.json</code>, <code className="font-mono text-[#dfba73] font-bold bg-[#090a0f] border border-[#222430] px-1 rounded">main.js</code>, and <code className="font-mono text-[#dfba73] font-bold bg-[#090a0f] border border-[#222430] px-1 rounded">README.md</code>) in a single directory named <code className="font-mono bg-[#090a0f] border border-[#222430] px-1 rounded text-white">elsepa-desktop/</code> on your computer.
            </li>
            <li>
              <strong>Install Dependencies</strong>: Open your terminal inside that directory and execute:
              <pre className="mt-2 bg-[#090a0f] text-[#33ff33]/85 border border-[#222430] p-2.5 rounded-lg font-mono text-[10px] select-all leading-normal">
                npm install
              </pre>
            </li>
            <li>
              <strong>Compile Portable Executable</strong>: Trigger the final bundler script:
              <div className="flex flex-col gap-1.5 mt-2">
                <span className="text-[10px] text-gray-400 font-bold">FOR WINDOWS (.exe):</span>
                <pre className="bg-[#090a0f] text-[#33ff33]/85 border border-[#222430] p-2 rounded-lg font-mono text-[10px] select-all leading-none">
                  npm run package-win
                </pre>
                <span className="text-[10px] text-gray-400 font-bold mt-1">FOR MACOS (.app):</span>
                <pre className="bg-[#090a0f] text-[#33ff33]/85 border border-[#222430] p-2 rounded-lg font-mono text-[10px] select-all leading-none">
                  npm run package-mac
                </pre>
              </div>
            </li>
          </ol>

          <div className="flex items-start gap-1 bg-amber-950/20 text-amber-300 border border-amber-900/30 p-3 rounded-lg text-xs leading-normal mt-2">
            <HelpCircle className="w-4 h-4 shrink-0 mt-0.5" />
            <div>
              <strong>Note:</strong> The downloaded configuration automatically sets the source URL to: <br />
              <span className="font-mono font-bold select-all bg-[#090a0f] border border-[#222430] px-1 rounded text-[10.5px] text-[#dfba73]">
                {template?.appUrl || "http://localhost:3000"}
              </span>
              . Any updates made in this workspace will sync directly with your desktop executable client!
            </div>
          </div>
        </div>

        {/* Tabulator code files preview */}
        <div className="lg:col-span-7 flex flex-col border border-[#222430] rounded-xl overflow-hidden h-[340px]">
          <div className="flex border-b border-[#222430] bg-[#090a0f]/45 p-1.5 select-none shrink-0">
            {["README.md", "package.json", "main.js"].map((f) => (
              <button
                key={f}
                id={`tab-file-${f}`}
                onClick={() => setActiveFile(f as any)}
                className={`px-3 py-1 text-xs font-bold rounded-lg cursor-pointer flex items-center gap-1.5 transition ${
                  activeFile === f
                    ? "bg-[#1c1d26] text-[#c5a059] border border-[#c5a059]/20 shadow-md"
                    : "text-[#94a3b8] hover:text-white"
                }`}
              >
                <FileCode className="w-3.5 h-3.5 text-gray-400" />
                {f}
              </button>
            ))}
          </div>

          <div className="flex-1 overflow-auto bg-[#090a0f]/95 p-4 font-mono text-[11px] text-gray-300 leading-normal scrollbar-thin">
            <pre className="whitespace-pre">
              {template ? template.files[activeFile] : "// Generating packaging templates..."}
            </pre>
          </div>
        </div>
      </div>
    </div>
  );
}
