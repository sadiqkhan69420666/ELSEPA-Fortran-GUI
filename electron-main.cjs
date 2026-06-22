const { app, BrowserWindow, Menu } = require('electron');
const path = require('path');
const { spawn } = require('child_process');
const http = require('http');

let mainWindow;
let serverProcess;

function startLocalServer() {
  const isWin = process.platform === 'win32';
  // Use the compiled production bundle if it exists, otherwise fallback to tsx dev mode
  const productionServerPath = path.join(__dirname, 'dist', 'server.cjs');
  const serverExists = require('fs').existsSync(productionServerPath);

  if (serverExists) {
    console.log('[Electron] Starting local production server from dist/server.cjs...');
    serverProcess = spawn('node', [productionServerPath], {
      cwd: __dirname,
      env: { ...process.env, NODE_ENV: 'production', PORT: '3000' }
    });
  } else {
    console.log('[Electron] Compiled bundle not found. Starting local development server...');
    const npxCmd = isWin ? 'npx.cmd' : 'npx';
    serverProcess = spawn(npxCmd, ['tsx', 'server.ts'], {
      cwd: __dirname,
      env: { ...process.env, NODE_ENV: 'development', PORT: '3000' }
    });
  }

  serverProcess.stdout.on('data', (data) => {
    console.log('[Server Standard Output]: ' + data);
  });

  serverProcess.stderr.on('data', (data) => {
    console.error('[Server Error Output]: ' + data);
  });
}

function pollLocalServerAndLoad(url, attempts = 0) {
  if (attempts > 40) {
    console.error('Failed to connect to local ELSEPA backend server on port 3000.');
    app.quit();
    return;
  }
  
  http.get(url, (res) => {
    if (res.statusCode === 200) {
      mainWindow.loadURL(url);
    } else {
      setTimeout(() => pollLocalServerAndLoad(url, attempts + 1), 250);
    }
  }).on('error', () => {
    setTimeout(() => pollLocalServerAndLoad(url, attempts + 1), 250);
  });
}

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1250,
    height: 850,
    title: "ELSEPA Desktop Simulation Lab",
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
    },
    backgroundColor: '#090a0f'
  });

  // Start back-end server automatically
  startLocalServer();

  // Load integrated elegant, dark, local loading indicator
  mainWindow.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(`
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
  `));

  // Poll local server and transition once the local port becomes responsive
  pollLocalServerAndLoad("http://localhost:3000");

  // Keep menu controls simple and clean
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
