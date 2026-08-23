// // gdc-docker-backend/server.js
// const express = require('express');
// const cors = require('cors');
// const { exec } = require('child_process');

// const app = express();
// app.use(cors({ origin: 'http://localhost:3000' }));

// // Start Docker app when backend starts
// exec('docker ps -q --filter "publish=8080" | xargs -r docker stop', () => {
//   exec('docker run -d -p 8081:80 genoxa', (error, stdout, stderr) => {
//     if (error) {
//       console.error('Failed to start Docker app:', stderr);
//     } else {
//       console.log('Docker app started with container ID:', stdout.trim());
//     }
//   });
// });

// const PORT = 3001;
// app.listen(PORT, () => {
//   console.log(`Backend running at http://localhost:${PORT}`);
// });

// gdc-docker-backend/server.js



const express = require('express');
const cors = require('cors');
const { exec } = require('child_process');

const app = express();
app.use(cors({ origin: 'http://localhost:3000' }));

// Image is built separately, by hand, in the docker-genoxa folder.
// This only runs the already-built image, so it works from any directory.
exec('docker rm -f genoxa', () => {
  exec("docker run -d --name genoxa -p 8081:8081 genoxa", (err, out, errText) => {
    if (err) console.error('Run failed:', errText);
    else console.log('Genoxa container started:', out.trim());
  });
});

const PORT = 3001;
app.listen(PORT, () => console.log(`Backend running at http://localhost:${PORT}`));