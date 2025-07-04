import { createSlice, PayloadAction } from "@reduxjs/toolkit";
import CohortBuilderDefaultConfig from "./data/cohort_builder.json";
import { CoreState } from "../../reducers";

export interface CohortBuilderCategoryConfig {
  readonly label: string;
  readonly queryOptions: {
    readonly docType: string;
    readonly indexType: string;
  };
  readonly facets: ReadonlyArray<string>;
}

type CohortBuilderCategory =
  | "general"
  | "demographic"
  | "general_diagnosis"
  | "disease_status"
  | "disease_specific_classifications"
  | "treatment"
  | "exposure"
  | "other_clinical_attributes"
  | "biospecimen"
  | "available_data"
  | "custom";

const initialState: { customFacets: string[] } = {
  customFacets: [],
};

export interface CohortBuilderCategoryFacet {
  readonly facetName: string;
}

const slice = createSlice({
  name: "cohort/builderConfig",
  initialState,
  reducers: {
    addFilterToCohortBuilder: (
      state,
      action: PayloadAction<CohortBuilderCategoryFacet>,
    ) => {
      if (!state.customFacets.includes(action.payload.facetName)) {
        state.customFacets = [...state.customFacets, action.payload.facetName];
      }
    },
    removeFilterFromCohortBuilder: (
      state,
      action: PayloadAction<CohortBuilderCategoryFacet>,
    ) => {
      state.customFacets = state.customFacets.filter(
        (x) => x != action.payload.facetName,
      );
    },
    resetCohortBuilderToDefault: (state) => {
      state.customFacets = [];
    },
  },
});

export const cohortBuilderConfigReducer = slice.reducer;
export const {
  addFilterToCohortBuilder,
  removeFilterFromCohortBuilder,
  resetCohortBuilderToDefault,
} = slice.actions;

export const selectCohortBuilderConfig = (
  state: CoreState,
): Record<CohortBuilderCategory, CohortBuilderCategoryConfig> => ({
  ...CohortBuilderDefaultConfig.config,
  custom: {
    ...CohortBuilderDefaultConfig.config.custom,
    facets: state.cohort.builderConfig.customFacets,
  },
});

/**
 * returns an array of all the filters used in the current configuration.
 * @param state - current core state/store
 */
export const selectCohortBuilderConfigFilters = (
  state: CoreState,
): string[] => [
  ...Object.values(CohortBuilderDefaultConfig.config).reduce(
    (filters: string[], category) => {
      return [...filters, ...category.facets];
    },
    [] as string[],
  ),
  ...state.cohort.builderConfig.customFacets,
];

export const selectCohortBuilderConfigCategory = (
  state: CoreState,
  category: CohortBuilderCategory,
): CohortBuilderCategoryConfig =>
  category === "custom"
    ? {
        ...CohortBuilderDefaultConfig.config.custom,
        facets: state.cohort.builderConfig.customFacets,
      }
    : CohortBuilderDefaultConfig.config[category];
