import path from 'node:path';
import { fileURLToPath } from 'node:url';
import express from 'express';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const app = express();
const port = 8080;

// static data
app.use('/data', express.static(path.join(__dirname, 'data')));

// // bundled assets
// app.use('/assets', express.static(path.join(__dirname, 'dist', 'assets')));

// app.use((req, res) => {
//   res.sendFile(path.join(__dirname, 'dist', 'index.html'));
// });

// start server
app.listen(port, () => {
  console.log(`PletzerBrowse running on port ${port}`);
});
