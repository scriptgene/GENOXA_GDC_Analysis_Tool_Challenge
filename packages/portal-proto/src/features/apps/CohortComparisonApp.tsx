import React, { useContext, useState } from "react";
import { useIsDemoApp } from "@/hooks/useIsDemoApp";
import {
  useCoreSelector,
  selectCurrentCohortName,
  FilterSet,
  selectCurrentCohortFilters,
  Cohort,
  selectCohortByIdOrName,
  selectCurrentCohortId,
  selectAllCohorts,
} from "@gff/core";
import { SelectionScreenContext } from "@gff/portal-components";
import CohortComparison from "../cohortComparison/CohortComparison";
import AdditionalCohortSelection from "@/features/cohortComparison/AdditionalCohortSelection";
import { useDeepCompareEffect } from "use-deep-compare";
import { usePrevious } from "@mantine/hooks";

export const cohortComparisonDemo1: {
  filter: FilterSet;
  name: string;
  id: string;
} = {
  filter: {
    mode: "and",
    root: {
      "cases.project.project_id": {
        operator: "includes",
        field: "cases.project.project_id",
        operands: ["TCGA-LGG"],
      },
      "genes.gene_id": {
        field: "genes.gene_id",
        operator: "includes",
        operands: ["ENSG00000138413", "ENSG00000182054"],
      },
    },
  },
  name: "Low grade gliomas - IDH1 or IDH2 mutated",
  id: "demoCohort1Id",
};

export const cohortComparisonDemo2: {
  filter: FilterSet;
  name: string;
  id: string;
} = {
  filter: {
    mode: "and",
    root: {
      "cases.project.project_id": {
        operator: "includes",
        field: "cases.project.project_id",
        operands: ["TCGA-LGG"],
      },
      "genes.gene_id": {
        field: "genes.gene_id",
        operator: "excludeifany",
        operands: ["ENSG00000138413", "ENSG00000182054"],
      },
      "cases.available_variation_data": {
        field: "cases.available_variation_data",
        operator: "includes",
        operands: ["ssm"],
      },
    },
  },
  name: "Low grade gliomas - IDH1 and IDH2 not mutated",
  id: "demoCohort2Id",
};

const CohortComparisonApp: React.FC = () => {
  const isDemoMode = useIsDemoApp();
  const { selectionScreenOpen, setSelectionScreenOpen, app, setActiveApp } =
    useContext(SelectionScreenContext);

  const allCohorts = useCoreSelector(selectAllCohorts);
  const allCohortsIds = Object.keys(allCohorts);

  /* Primary Cohort Details */
  const primaryCohortName = useCoreSelector((state) =>
    selectCurrentCohortName(state),
  );
  const primaryCohortId = useCoreSelector((state) =>
    selectCurrentCohortId(state),
  );
  const primaryCohortFilter = useCoreSelector((state) =>
    selectCurrentCohortFilters(state),
  );
  /* Primary Cohort Details End */

  /* Comparison Cohort Details */
  const [comparisonCohort, setComparisonCohort] = useState<Cohort>();
  const comparisonCohortId = comparisonCohort?.id;
  const comparisonCohortObj: Cohort = useCoreSelector((state) =>
    selectCohortByIdOrName(state, comparisonCohort?.id, comparisonCohort?.name),
  );
  const comparisonCohortFilter = comparisonCohortObj?.filters;
  /* Comparison Cohort Details End */

  const cohorts = isDemoMode
    ? {
        primary_cohort: cohortComparisonDemo1,
        comparison_cohort: cohortComparisonDemo2,
      }
    : {
        primary_cohort: {
          filter: primaryCohortFilter,
          name: primaryCohortName,
          id: primaryCohortId,
        },
        comparison_cohort: {
          filter: comparisonCohortFilter,
          name: comparisonCohort?.name,
          id: comparisonCohort?.id,
        },
      };

  const prevPrimaryCohortId = usePrevious(primaryCohortId);
  const prevComparisonCohortId = usePrevious(comparisonCohortId);

  useDeepCompareEffect(() => {
    if (
      !isDemoMode &&
      (!allCohortsIds.includes(prevPrimaryCohortId) ||
        !allCohortsIds.includes(prevComparisonCohortId))
    ) {
      setSelectionScreenOpen(true);
    }
  }, [
    isDemoMode,
    allCohortsIds,
    prevPrimaryCohortId,
    prevComparisonCohortId,
    setSelectionScreenOpen,
  ]);

  return selectionScreenOpen ? (
    <AdditionalCohortSelection
      app={app}
      setOpen={setSelectionScreenOpen}
      setActiveApp={setActiveApp}
      setComparisonCohort={setComparisonCohort}
    />
  ) : (
    <CohortComparison cohorts={cohorts} demoMode={isDemoMode} />
  );
};

export default CohortComparisonApp;
