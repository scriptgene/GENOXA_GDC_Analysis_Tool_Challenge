// gdc-docker-backend/server.js
const express = require('express');
const cors = require('cors');
const { exec } = require('child_process');

const app = express();
app.use(cors({ origin: 'http://localhost:3000' }));

// Start Docker app when backend starts
exec('docker ps -q --filter "publish=8080" | xargs -r docker stop', () => {
  exec('docker run -d -p 8081:80 genoxa', (error, stdout, stderr) => {
    if (error) {
      console.error('Failed to start Docker app:', stderr);
    } else {
      console.log('Docker app started with container ID:', stdout.trim());
    }
  });
});

const PORT = 3001;
app.listen(PORT, () => {
  console.log(`Backend running at http://localhost:${PORT}`);
});