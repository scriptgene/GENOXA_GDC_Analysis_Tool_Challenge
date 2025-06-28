# GENOXA: A Tool for Analyzing Bulk RNA-Seq Data and Building Machine Learning Models with Differentially Expressed Genes

### TCGA Data from GDC Data Portal Integrated with GENOXA Tool

Bulk RNA-Seq data for the TCGA project were downloaded from the [GDC data portal](https://portal.gdc.cancer.gov/).  
The dataset falls under the **Transcriptome Profiling** category, utilizing the **RNA-Seq** experimental strategy and the **STAR - Counts** workflow. 
It includes **Primary Tumor** and **Solid Tissue Normal** sample types and is publicly accessible with **open access**.

Each dataset includes **molecular, clinical, and expression data** for different cancer types, allowing for in-depth research and analysis.

#### Available Datasets

1. **TCGA-PRAD** *(Prostate Adenocarcinoma)* – Prostate Cancer  
2. **TCGA-BRCA** *(Breast Invasive Carcinoma)* – Breast Cancer  
3. **TCGA-LUAD** *(Lung Adenocarcinoma)* – Lung Cancer  
4. **TCGA-KIRC** *(Kidney Renal Clear Cell Carcinoma)* – Kidney Cancer  
5. **TCGA-COAD** *(Colon Adenocarcinoma)* – Colon Cancer  
6. **TCGA-STAD** *(Stomach Adenocarcinoma)* – Stomach Cancer  
7. **TCGA-BLCA** *(Bladder Urothelial Carcinoma)* – Bladder Cancer  
8. **TCGA-LIHC** *(Liver Hepatocellular Carcinoma)* – Liver Cancer

## 1. GDC Frontend Framework Installation

### Prerequisites

This is a multi-workspace repo that requires npm v10.2.4. The minimum node version is set to v20.11.0.

Node can be downloaded from the official Node.js site. You may also consider using a [Node version manager](https://docs.npmjs.com/cli/v7/configuring-npm/install#using-a-node-version-manager-to-install-nodejs-and-npm).

Your version of Node may not ship with npm v10.2.4. To install it, run:

```bash
npm install npm@10.2.4
```

If you are using a Node version manager, you can run the following to install the correct version of Node:

```bash
nvm install 20.11.0
```

to use the correct version of Node:

```bash
nvm use 20.11.0
```

### Install Dependencies

From the root of the project, install dependencies by running:

```bash
npm install
```

Installing from the root of the repository is required to avoid
multiple installations of React in the workspaces. When this happens,
React will fail to render.

### Adding Dependencies

Dependencies can also be installed from the root of the repository.
Error might occured for magin string module, run this command:-

```bash
npm install @napi-rs/magic-string-darwin-x64
```

Since this is a TypeScript project, installing the community type definitions may also be required:

```bash
npm install --save-dev @types/my-package --workspace=packages/core
```

## 2. Create docker image of GENOXA App

### Install Docker

Before running this project, make sure you have Docker installed on your system.

Build docker image by going into the specific directory --> gc-docker-backend --> NCI_GENOXA_Lite
```bash
docker build -t genoxa .
```

Run the backend server by going into root --> gdc-docker-backend, then build the docker image and then run the docker image
```bash
node server.js
```
## 3. Development

Run the prototype in dev mode with auto-rebuilding into new terminal:

```bash
npm run dev
```

By default, this will start a dev server listening to http://localhost:3000

## 🖼️ Screenshot

![App Screenshot](./assets/1.png)
![App Screenshot](./assets/2.png)
![App Screenshot](./assets/3.png)
![App Screenshot](./assets/4.png)
![App Screenshot](./assets/5.png)
![App Screenshot](./assets/6.png)


## 🎬 Video tutorial 

[Click here to watch the video of GENOXA App integrated with GDC data portal](https://youtu.be/WNrVuYOJjtI)

## 🎬 GENOXA App only

[Click here to watch the video of GENOXA App only](https://www.youtube.com/watch?v=6nPyV5qtOUc)

