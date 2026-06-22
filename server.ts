import express from "express";
import path from "path";
import cors from "cors";
import { createServer as createViteServer } from "vite";
import dotenv from "dotenv";
import { execSync } from "child_process";
import fs from "fs";
import AdmZip from "adm-zip";
import { runScatteringSimulation, PERIODIC_TABLE } from "./src/utils/physicsSolver";
import { SimulationParams, DataPoint, PhaseShiftPoint } from "./src/types";

// Load environment variables
dotenv.config();

// Model code maps for Fortran program arguments
const potentialModelMap = {
  "dirac-fock": 1,
  "hartree-fock": 2,
  "bohr-screening": 3,
  "yukawa": 4
};

const nuclearModelMap = {
  "point": 1,
  "uniform": 2,
  "fermi": 3
};

const exchangeModelMap = {
  "none": 0,
  "furness-mccarthy": 1,
  "riley-truhlar": 2
};

function parseElsepaOut(filePath: string): {
  deBroglieWavelength: number,
  totalElasticCrossSection: number,
  phaseShifts: PhaseShiftPoint[],
  dcsData: DataPoint[]
} | null {
  try {
    if (!fs.existsSync(filePath)) return null;
    const content = fs.readFileSync(filePath, "utf-8");
    const lines = content.split("\n");

    let deBroglieWavelength = 0;
    let totalElasticCrossSection = 0;
    const phaseShifts: PhaseShiftPoint[] = [];
    const dcsData: DataPoint[] = [];

    let readingPhaseShifts = false;
    let readingAngular = false;

    for (const rawLine of lines) {
      const line = rawLine.trim();
      if (!line) continue;

      if (line.startsWith("de Broglie Wavelength (A):")) {
        deBroglieWavelength = parseFloat(line.split(":")[1].trim());
      } else if (line.startsWith("Elastic Cross Section (a0^2):")) {
        totalElasticCrossSection = parseFloat(line.split(":")[1].trim());
      } else if (line.startsWith("PHASE SHIFTS:")) {
        readingPhaseShifts = true;
        readingAngular = false;
        continue;
      } else if (line.startsWith("ANGULAR DISTRIBUTIONS:")) {
        readingPhaseShifts = false;
        readingAngular = true;
        continue;
      }

      // Check column header skips
      if (line.startsWith("l ") || line.startsWith("Angle(deg)") || line.startsWith("l_") || line.startsWith("====")) {
        continue;
      }

      if (readingPhaseShifts) {
        // Stop if we hit any footer lines
        if (line.includes("====") || line.startsWith("ANGULAR")) {
          readingPhaseShifts = false;
        } else {
          const parts = line.split(/\s+/);
          if (parts.length >= 3) {
            const lVal = parseInt(parts[0]);
            const deltaVal = parseFloat(parts[1]);
            const etaVal = parseFloat(parts[2]);
            if (!isNaN(lVal)) {
              phaseShifts.push({ l: lVal, delta: deltaVal, eta: etaVal });
            }
          }
        }
      }

      if (readingAngular) {
        if (line.includes("====")) {
          readingAngular = false;
        } else {
          const parts = line.split(/\s+/);
          if (parts.length >= 4) {
            const angleVal = parseFloat(parts[0]);
            const dcsElastic = parseFloat(parts[1]);
            const dcsRuth = parseFloat(parts[2]);
            const sherman = parseFloat(parts[3]);
            if (!isNaN(angleVal)) {
              dcsData.push({
                angle: angleVal,
                dcs: dcsElastic,
                dcsRutherford: dcsRuth,
                Sherman: sherman
              });
            }
          }
        }
      }
    }

    if (dcsData.length === 0) return null;

    return {
      deBroglieWavelength,
      totalElasticCrossSection,
      phaseShifts,
      dcsData
    };
  } catch (err) {
    console.error("Error parsing elsepa.out file:", err);
    return null;
  }
}

async function startServer() {
  const app = express();
  const PORT = 3000;

  // Body parsing and security middleware
  app.use(express.json({ limit: "15mb" }));
  app.use(cors());

  // API Route: Health check
  app.get("/api/health", (req, res) => {
    res.json({ status: "ok", time: new Date().toISOString() });
  });

  // API Route: Run Scattering Simulation (Local Fortran with TS Fallback)
  app.post("/api/simulate", (req, res) => {
    try {
      const params = req.body as SimulationParams;
      if (!params || !params.atomicNumber || !params.energy) {
        return res.status(400).json({ error: "Missing required simulation parameters." });
      }

      const element = PERIODIC_TABLE.find(el => el.number === params.atomicNumber) || PERIODIC_TABLE[5];
      const energyInEv = params.energyUnit === "keV" ? params.energy * 1000 : params.energyUnit === "MeV" ? params.energy * 1000000 : params.energy;

      // 1. Try to see if gfortran compiler is available
      let gfortranAvailable = false;
      try {
        execSync("gfortran --version", { stdio: "ignore" });
        gfortranAvailable = true;
      } catch {
        // no gfortran
      }

      if (gfortranAvailable) {
        console.log("[Fortran Engine] Found gfortran compiler. Attempting native run.");
        // Generate elsepa.in
        const inputContent = [
          params.atomicNumber,
          params.projectile === "electron" ? 1 : -1,
          energyInEv,
          potentialModelMap[params.potentialModel] || 1,
          nuclearModelMap[params.nuclearModel] || 1,
          exchangeModelMap[params.exchangeModel] || 0,
          params.absorptionModel ? 1 : 0,
          params.correlationPolarization ? 1 : 0
        ].join("\n") + "\n";

        fs.writeFileSync("elsepa.in", inputContent);

        // Confirm compiled binary exists, compile if needed
        const isWin = process.platform === "win32";
        const binaryName = isWin ? "elsepa_solver.exe" : "./elsepa_solver";
        const binaryPath = path.join(process.cwd(), isWin ? "elsepa_solver.exe" : "elsepa_solver");

        if (!fs.existsSync(binaryPath)) {
          console.log("[Fortran Engine] Compiling elsepa_solver.f90...");
          const compileCmd = isWin
            ? "gfortran elsepa_solver.f90 -o elsepa_solver.exe"
            : "gfortran elsepa_solver.f90 -o elsepa_solver";
          execSync(compileCmd, { cwd: process.cwd(), stdio: "inherit" });
        }

        // Run Fortran binary
        execSync(binaryName, { cwd: process.cwd(), timeout: 15000 });

        // Parse elsepa.out
        const parsed = parseElsepaOut("elsepa.out");
        if (parsed) {
          const elsepaInFile = fs.existsSync("elsepa.in") ? fs.readFileSync("elsepa.in", "utf-8") : "";
          const elsepaOutFile = fs.existsSync("elsepa.out") ? fs.readFileSync("elsepa.out", "utf-8") : "";

          console.log("[Fortran Engine] Real Fortran simulation successfully parsed!");
          return res.json({
            params,
            element,
            dcsData: parsed.dcsData,
            phaseShifts: parsed.phaseShifts,
            totalElasticCrossSection: parsed.totalElasticCrossSection,
            momentumTransferCrossSection: parsed.totalElasticCrossSection * 0.9, // placeholder momentum integration
            deBroglieWavelength: parsed.deBroglieWavelength,
            elsepaInFile,
            elsepaOutFile,
            timestamp: new Date().toISOString(),
            engine: "native_fortran"
          });
        }
      }

      // 2. Gracious Fallback to TS Mathematical simulation if Fortran tool-chain is absent
      console.log("[Fortran Engine] gfortran not found. Falling back to high-fidelity JS/TS emulator.");
      const result = runScatteringSimulation(params);
      return res.json({
        ...result,
        engine: "ts_emulated"
      });

    } catch (err: any) {
      console.error("[Fortran Engine Error] Failed run:", err);
      // Fallback in case of code or file system runtime crash
      try {
        const result = runScatteringSimulation(req.body);
        return res.json({
          ...result,
          engine: "ts_emulated",
          errorMsg: err.message
        });
      } catch (fallbackErr: any) {
        return res.status(500).json({ error: err.message, fallbackError: fallbackErr.message });
      }
    }
  });

  // API Route: Serve pre-packaged Desktop Wrapper configuration (JSON/scripts)
  app.get("/api/desktop/template", (req, res) => {
    // Serve a JSON configuration with all files needed for the desktop wrapper that runs completely locally
    res.json({
      appName: "ELSEPA Desktop Lab",
      appUrl: "http://localhost:3000",
      files: {
        "package.json": JSON.stringify({
          name: "elsepa-desktop-client",
          version: "1.0.0",
          description: "Portable Independent Local Desktop Client for ELSEPA Simulation Suite",
          main: "main.js",
          scripts: {
            "start": "electron .",
            "package-win": "npx electron-packager . ELSEPA-Suite --platform=win32 --arch=x64 --out=dist --overwrite --icon=icon",
            "package-mac": "npx electron-packager . ELSEPA-Suite --platform=darwin --arch=x64 --out=dist --overwrite --icon=icon",
            "package-linux": "npx electron-packager . ELSEPA-Suite --platform=linux --arch=x64 --out=dist --overwrite"
          },
          dependencies: {
            "electron": "^28.2.0"
          },
          devDependencies: {
            "electron-packager": "^17.1.2"
          }
        }, null, 2),
        "main.js": `const { app, BrowserWindow, Menu } = require('electron');
const path = require('path');
const { spawn } = require('child_process');
const http = require('http');

let mainWindow;
let serverProcess;

function startLocalServer() {
  const isWin = process.platform === 'win32';
  const fs = require('fs');
  const serverPath = path.join(__dirname, 'dist', 'server.cjs');
  const serverExists = fs.existsSync(serverPath);

  if (serverExists) {
    console.log('[Electron] Starting local production server from dist/server.cjs...');
    serverProcess = spawn('node', [serverPath], {
      cwd: __dirname,
      env: { ...process.env, NODE_ENV: 'production', PORT: '3000' }
    });
  } else {
    console.log('[Electron] Compiled production bundle not found. Starting local development server using tsx...');
    const npxCmd = isWin ? 'npx.cmd' : 'npx';
    serverProcess = spawn(npxCmd, ['tsx', 'server.ts'], {
      cwd: __dirname,
      env: { ...process.env, NODE_ENV: 'development', PORT: '3000' }
    });
  }

  serverProcess.stdout.on('data', (data) => {
    console.log('[Local Server Output]: ' + data);
  });

  serverProcess.stderr.on('data', (data) => {
    console.error('[Local Server Error]: ' + data);
  });
}

function pollLocalServerAndLoad(url, attempts = 0) {
  if (attempts > 60) {
    console.error('Failed to connect to local ELSEPA backend server on port 3000.');
    app.quit();
    return;
  }
  
  http.get(url, (res) => {
    if (res.statusCode === 200) {
      mainWindow.loadURL(url);
    } else {
      setTimeout(() => pollLocalServerAndLoad(url, attempts + 1), 200);
    }
  }).on('error', () => {
    setTimeout(() => pollLocalServerAndLoad(url, attempts + 1), 200);
  });
}

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1280,
    height: 850,
    title: "ELSEPA Desktop Simulation Lab",
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      sandbox: true
    },
    backgroundColor: '#090a0f'
  });

  // Start the offline local backend automatically
  startLocalServer();

  // Load integrated elegant offline loading indicator
  mainWindow.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(\`
    <!DOCTYPE html>
    <html>
    <head>
      <title>ELSEPA Lab Loading</title>
      <style>
        body {
          background-color: #090a0f;
          color: #94a3b8;
          font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          height: 100vh;
          margin: 0;
        }
        .spinner {
          border: 3px solid #222430;
          border-top: 3px solid #c5a059;
          border-radius: 50%;
          width: 36px;
          height: 36px;
          animation: spin 1s linear infinite;
          margin-bottom: 20px;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        h2 {
          color: white;
          font-weight: 500;
          margin-bottom: 8px;
        }
        p {
          font-size: 13px;
        }
      </style>
    </head>
    <body>
      <div class="spinner"></div>
      <h2>Starting Native ELSEPA Desktop Lab</h2>
      <p>Initializing offline-independent physics solver server on local port 3000...</p>
    </body>
    </html>
  \`));

  // Poll local server and present the full UI inside Electron once port is ready with no sign-ins!
  pollLocalServerAndLoad("http://localhost:3000");

  Menu.setApplicationMenu(null);

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

app.whenReady().then(createWindow);

app.on('window-all-closed', () => {
  if (serverProcess) {
    serverProcess.kill();
  }
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});
`,
        "README.md": `# ELSEPA Desktop Suite - Offline Independent Wrapper

This directory contains the necessary Electron configuration to package the **ELSEPA Online Simulation Suite** into a completely local, offline-independent desktop application for Windows, macOS, or Linux.

This configuration automatically launches the local Node.js Express server on port 3000 upon startup, compiles and runs the native Fortran binaries, and renders the frontend with zero connections to the cloud. **There are no Google accounts, sign-ups, or internet dependencies required.**

## Requirements
- Node.js (with npm) installed on your system.
- \`gfortran\` compiler (recommended, to run native Dirac/Schrödinger solver instead of TS emulator).

## Quick Setup & Start

1. **Install dependencies**:
   \`\`\`bash
   npm install
   \`\`\`

2. **Launch the offline client**:
   \`\`\`bash
   npm start
   \`\`\`

3. **Package as portable executable files**:
   - **For Windows**:
     \`\`\`bash
     npm run package-win
     \`\`\`
   - **For macOS**:
     \`\`\`bash
     npm run package-mac
     \`\`\`
   - **For Linux**:
     \`\`\`bash
     npm run package-linux
     \`\`\`

Your executable is saved under \`dist/\`. Copy it anywhere and execute it instantly with zero friction!
`
      }
    });
  });

  // API Route: Download complete source directory as ZIP for fully-offline local Electron lab
  app.get("/api/desktop/download-zip", (req, res) => {
    try {
      const zip = new AdmZip();
      const rootDir = process.cwd();

      // Add single files
      const singleFiles = [
        "package.json",
        "vite.config.ts",
        "tsconfig.json",
        "index.html",
        "server.ts",
        "electron-main.cjs",
        ".gitignore",
        "elsepa_solver.f90"
      ];

      for (const file of singleFiles) {
        const filePath = path.join(rootDir, file);
        if (fs.existsSync(filePath)) {
          zip.addLocalFile(filePath);
        }
      }

      // Add folders recursively
      const folders = ["src", "assets"];
      for (const folder of folders) {
        const folderPath = path.join(rootDir, folder);
        if (fs.existsSync(folderPath)) {
          zip.addLocalFolder(folderPath, folder);
        }
      }

      const zipBuffer = zip.toBuffer();
      res.set({
        "Content-Type": "application/zip",
        "Content-Disposition": "attachment; filename=elsepa-desktop-lab.zip",
        "Content-Length": zipBuffer.length
      });
      res.send(zipBuffer);
    } catch (err: any) {
      console.error("Failed to generate offline zip:", err);
      res.status(500).json({ error: "Failed to generate ZIP archive: " + err.message });
    }
  });

  // Serve static files in production, use Vite middleware in dev
  if (process.env.NODE_ENV !== "production") {
    const vite = await createViteServer({
      server: { middlewareMode: true },
      appType: "spa",
    });
    app.use(vite.middlewares);
  } else {
    const distPath = path.join(process.cwd(), "dist");
    app.use(express.static(distPath));
    app.get("*", (req, res) => {
      res.sendFile(path.join(distPath, "index.html"));
    });
  }

  app.listen(PORT, "0.0.0.0", () => {
    console.log(`[ELSEPA Suite Server] running on http://localhost:${PORT}`);
  });
}

startServer().catch((err) => {
  console.error("Failed to start server:", err);
  process.exit(1);
});
