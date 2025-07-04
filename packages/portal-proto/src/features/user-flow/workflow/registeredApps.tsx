import SelectionPanel from "@/features/set-operations/SelectionPanel";
import AdditionalCohortSelection from "@/features/cohortComparison/AdditionalCohortSelection";
import DownloadAllButton from "@/features/cDave/DownloadAllButton";
import SearchInput from "@/components/SearchInput";
import ClinicalDataIcon from "public/user-flow/icons/apps/ClinicalData.svg";
import CohortBuilderIcon from "public/user-flow/icons/apps/CohortBuilder.svg";
import MutationFrequencyIcon from "public/user-flow/icons/apps/MutationsFrequency.svg";
import RepositoryIcon from "public/user-flow/icons/database.svg";
import ProjectsIcon from "public/user-flow/icons/crowd-of-users.svg";
import CohortComparisonIcon from "public/user-flow/icons/apps/CohortComparison.svg";
import SetOperationsIcon from "public/user-flow/icons/apps/SetOperations.svg";
import SequenceReadsIcon from "public/user-flow/icons/apps/SequenceReads.svg";
import BAMSlicingDownloadIcon from "public/user-flow/icons/apps/BAMSlicingDownload.svg";
import CohortLevelMAFIcon from "public/user-flow/icons/apps/CohortLevelMAF.svg";
import ProteinPaintIcon from "public/user-flow/icons/apps/ProteinPaint.svg";
import OncoMatrixIcon from "public/user-flow/icons/apps/OncoMatrix.svg";
import GeneExpressionIcon from "public/user-flow/icons/apps/GeneExpression.svg";
import ScRNASeqIcon from "public/user-flow/icons/apps/scRNASeq.svg";
import CopyNumberSegmentIcon from "public/user-flow/icons/apps/CopyNumberSegment.svg";
import { useLazyScRNAseqCaseCountQuery } from "../../proteinpaint/scRNAseqCaseCount";
import { useLazyCnvSegmentCaseCountQuery } from "../../proteinpaint/CnvSegmentCaseCount";
import { CountHookRegistry } from "@gff/core";
import { AppRegistrationEntry } from "@gff/portal-components";

CountHookRegistry.getInstance().registerHook(
  "scRNAseqCaseCount",
  useLazyScRNAseqCaseCountQuery,
);

CountHookRegistry.getInstance().registerHook(
  "CnvSegmentCaseCount",
  useLazyCnvSegmentCaseCountQuery,
);

export const COHORTS = [
  { name: "New Custom Cohort", facets: [] },
  {
    name: "Baily's Cohort",
    facets: [
      { field: "cases.primary_site", value: ["breast", "bronchus and lung"] },
    ],
  },
  {
    name: "Pancreas",
    facets: [{ field: "cases.primary_site", value: ["pancreas"] }],
  },
];

export const REGISTERED_APPS: AppRegistrationEntry[] = [
  {
    name: "Clinical Data Analysis",
    icon: <ClinicalDataIcon aria-hidden="true" />,
    tags: ["clinicalAnalysis"],
    hasDemo: true,
    countsField: "caseCount",
    id: "CDave",
    description:
      "Use clinical variables to perform basic statistical analysis of your cohort.",
    noDataTooltip:
      "Current cohort does not have cases available for visualization.",
    rightComponent: DownloadAllButton,
  },
  {
    name: "Cohort Builder",
    icon: (
      <CohortBuilderIcon
        width={64}
        height={64}
        viewBox="0 0 60 60"
        aria-hidden="true"
      />
    ),
    tags: ["generalUtility"],
    hasDemo: false,
    id: "CohortBuilder",
    countsField: "caseCount",
    description:
      "Build and define your custom cohorts using a variety of clinical and biospecimen features.",
    rightComponent: SearchInput,
  },
  {
    name: "Mutation Frequency",
    icon: <MutationFrequencyIcon />,
    tags: ["variantAnalysis", "ssm"],
    hasDemo: true,
    id: "MutationFrequencyApp",
    countsField: "cnvOrSsmCaseCount",
    description:
      "Visualize most frequently mutated genes and somatic mutations.",
    noDataTooltip:
      "Current cohort does not have SSM or CNV data available for visualization.",
    optimizeRules: ["something == something"],
  },
  {
    name: "Repository",
    icon: <RepositoryIcon aria-hidden="true" />,
    tags: ["files"],
    hasDemo: false,
    countsField: "caseCount",
    id: "Downloads",
    description:
      "Browse and download the files associated with your cohort for more sophisticated analysis.",
  },
  {
    name: "Projects",
    icon: (
      <ProjectsIcon
        width={64}
        height={64}
        viewBox="0 -20 128 128"
        aria-hidden="true"
      />
    ),
    tags: [],
    hasDemo: false,
    id: "Projects",
    countsField: "caseCount",
    description:
      "View the Projects available within the GDC and select them for further exploration and analysis.",
  },
  {
    name: "Cohort Comparison",
    icon: <CohortComparisonIcon aria-hidden="true" />,
    tags: ["clinicalAnalysis"],
    hasDemo: true,
    id: "CohortComparisonApp",
    countsField: "caseCount",
    description:
      "Display the survival analysis of your cohorts and compare characteristics such as gender, vital status and age at diagnosis.",
    noDataTooltip:
      "Current cohort does not have cases available for visualization.",
    selectionScreen: AdditionalCohortSelection,
  },
  // TODO uncomment to show gene expression card
  /*
  {
    name: "Gene Expression",
    icon: "icons/aops/GeneExpression.svg",
    tags: ["geneExpression"],
    hasDemo: true,
    id: "GeneExpression",
    description: "Visualize patterns in gene expression in your cohort.",
    caseCounts: 0.11,
    optimizeRules: ["something == something"],
    noDataTooltip:
      "Current cohort does not have gene expression data available for visualization.",
  },
  */
  {
    name: "Set Operations",
    icon: <SetOperationsIcon aria-hidden="true" />,
    tags: ["generalUtility"],
    hasDemo: true,
    hideCounts: true,
    countsField: "caseCount",
    description:
      "Display a Venn diagram and compare/contrast your cohorts or sets of the same type.",
    id: "SetOperationsApp",
    selectionScreen: SelectionPanel,
  },
  {
    name: "Sequence Reads",
    icon: <SequenceReadsIcon aria-hidden="true" />,
    tags: ["sequenceAnalysis"],
    hasDemo: false,
    countsField: "sequenceReadCaseCount",
    description:
      "Visualize sequencing reads for a given gene, position, SNP, or variant.",
    id: "SequenceReadApp",
    noDataTooltip:
      "Current cohort does not have available BAMs for visualization.",
    optimizeRules: ["data format = BAM"],
  },
  {
    name: "BAM Slicing Download",
    icon: <BAMSlicingDownloadIcon aria-hidden="true" />,
    tags: ["sequenceAnalysis"],
    hasDemo: false,
    countsField: "sequenceReadCaseCount",
    description: "Download a BAM slice.",
    id: "BamDownloadApp",
    noDataTooltip: "Current cohort does not have available BAMs for download.",
    optimizeRules: ["data format = BAM"],
  },
  {
    name: "ProteinPaint",
    icon: (
      <ProteinPaintIcon
        height={48}
        width={80}
        viewBox="-12 0 80 48"
        aria-hidden="true"
      />
    ),
    tags: ["variantAnalysis", "ssm"],
    hasDemo: true,
    description:
      "Visualize mutations in protein-coding genes by consequence type and protein domain.",
    id: "ProteinPaintApp",
    countsField: "ssmCaseCount",
    optimizeRules: ["available data = ssm"],
    noDataTooltip:
      "Current cohort does not have SSM data available for visualization.",
  },
  {
    name: "OncoMatrix",
    icon: (
      <OncoMatrixIcon
        className="m-auto"
        height={48}
        width={80}
        aria-hidden="true"
      />
    ),
    tags: ["variantAnalysis", "cnv", "ssm"],
    hasDemo: true,
    description:
      "Visualize the top most mutated cases and genes affected by high impact mutations in your cohort.",
    id: "OncoMatrix",
    countsField: "cnvOrSsmCaseCount",
    optimizeRules: ["available data = ssm or cnv"],
    noDataTooltip:
      "Current cohort does not have SSM or CNV data available for visualization.",
  },
  {
    name: "Single Cell RNA-seq",
    icon: (
      <ScRNASeqIcon
        className="m-auto"
        height={48}
        width={80}
        aria-hidden="true"
      />
    ),
    tags: ["variantAnalysis", "cnv", "ssm"],
    hasDemo: true,
    description:
      "Visualize single-cell RNA-Seq data with cluster plots and gene expression overlays.",
    id: "scRNAseq",
    countsField: "scRNAseqCaseCount",
    optimizeRules: ["available data = ssm or cnv"],
    noDataTooltip:
      "Current cohort does not have scRNAseq data available for visualization.",
  },
  {
    name: "Gene Expression Clustering",
    icon: (
      <GeneExpressionIcon
        className="m-auto"
        height={48}
        width={80}
        aria-hidden="true"
      />
    ),
    tags: ["variantAnalysis", "cnv", "ssm"],
    hasDemo: true,
    description:
      "Visualize the top most variably expressed genes in your cohort.",
    id: "GeneExpression",
    countsField: "geneExpressionCaseCount",
    optimizeRules: ["available data = ssm or cnv"],
    noDataTooltip:
      "Current cohort does not have gene expression data available for visualization.",
  },
  {
    name: "Cohort Level MAF",
    // TODO: replace with a MAF icon
    icon: <CohortLevelMAFIcon aria-hidden="true" />,
    tags: ["ssm"],
    hasDemo: false,
    countsField: "mafFileCount",
    description: "Combine and download cohort level MAF data.",
    id: "CohortLevelMAF",
    noDataTooltip:
      "Current cohort does not have MAF data available for download.",
  },
  {
    name: "Copy Number Segment",
    icon: (
      <CopyNumberSegmentIcon
        height={48}
        width={80}
        viewBox="-12 0 80 48"
        aria-hidden="true"
      />
    ),
    tags: ["variantAnalysis", "cnv"],
    hasDemo: true,
    description: "Visualize copy number variation over a gene or region.",
    id: "CnvSegmentApp",
    countsField: "CnvSegmentCaseCount",
    noDataTooltip:
      "Current cohort does not have CNV segment data available for visualization.",
  },
  {
    name: "GENOXA",
    icon: <GeneExpressionIcon aria-hidden="true" />,
    tags: ["clinicalAnalysis"],
    hasDemo: false,
    countsField: "caseCount",
    hideCounts: true,
    id: "SagarPatel_GENOXA",
    description:
      "Visualize differentially expressed (DE) genes with machine learning (ML) based modeling",
    noDataTooltip:
      "Current cohort does not have cases available for visualization.",
    rightComponent: DownloadAllButton,
  },
];

export const APPTAGS = [
  { value: "clinicalAnalysis", name: "Clinical Analysis" },
  { value: "generalUtility", name: "General Utility" },
  { value: "variantAnalysis", name: "Variant Analysis" },
  // TODO uncomment to show gene expression tag
  // { value: "geneExpression", name: "Gene Expression" },
  { value: "sequenceAnalysis", name: "Sequence Analysis" },
  { value: "cnv", name: "CNV" },
  { value: "ssm", name: "SSM" },
];

export const RECOMMENDED_APPS = ["Projects", "CohortBuilder", "Downloads"];
