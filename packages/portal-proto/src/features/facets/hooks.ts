import { useEffect, useMemo } from "react";
import { UnknownAction, ThunkDispatch } from "@reduxjs/toolkit";
import { useDeepCompareEffect } from "use-deep-compare";
import {
  Operation,
  EnumValueExtractorHandler,
  EnumOperandValue,
  OperandValue,
  FacetBuckets,
  handleOperation,
  FilterSet,
  removeCohortFilter,
  fetchFacetByNameGQL,
  updateCohortFilter,
  useCoreDispatch,
  useCoreSelector,
  GQLDocType,
  GQLIndexType,
  NumericFromTo,
  selectRangeFacetByField,
  fetchFacetContinuousAggregation,
  selectFacetByDocTypeAndField,
  usePrevious,
  selectTotalCountsByName,
  selectMultipleFacetsByDocTypeAndField,
  selectCurrentCohortFiltersByName,
  selectCurrentCohortFiltersByNames,
  GqlOperation,
  buildCohortGqlOperator,
  selectCurrentCohortFilters,
  useCurrentCohortFilters,
  CoreState,
  selectCohortBuilderConfig,
  selectFacetDefinition,
} from "@gff/core";
import isEqual from "lodash/isEqual";
import { getFacetInfo } from "@/features/cohortBuilder/utils";

/**
 * Filter selector for all the facet filters
 */
const useCohortFacetFilter = (): FilterSet => {
  return useCoreSelector((state) => selectCurrentCohortFilters(state));
};

export const extractValue = (op: Operation): EnumOperandValue => {
  const handler = new EnumValueExtractorHandler();
  return handleOperation<EnumOperandValue>(handler, op);
};

/**
 * Selector for the facet values (if any) from the current cohort
 * @param field - field name to find filter for
 * @returns Value of Filters or undefined
 */
const useCohortFacetFilterByName = (field: string): OperandValue => {
  const enumFilters: Operation = useCoreSelector((state) =>
    selectCurrentCohortFiltersByName(state, field),
  );
  return enumFilters ? extractValue(enumFilters) : undefined;
};

const useEnumFiltersByNames = (
  fields: ReadonlyArray<string>,
): Record<string, OperandValue> => {
  const enumFilters: Record<string, Operation> = useCoreSelector((state) =>
    selectCurrentCohortFiltersByNames(state, fields),
  );
  return Object.entries(enumFilters).reduce((obj, [key, value]) => {
    if (value) obj[key] = extractValue(value);
    return obj;
  }, {});
};

/**
 *  Facet Selector using GQL which will refresh when filters/enum values changes.
 */
export const useEnumFacet = (
  field: string,
  { docType, indexType }: { docType: GQLDocType; indexType: GQLIndexType },
) => {
  const coreDispatch = useCoreDispatch();

  const facet: FacetBuckets = useCoreSelector((state) =>
    selectFacetByDocTypeAndField(state, docType, field),
  );

  const enumValues = useCohortFacetFilterByName(field);
  const currentCohortFilters = useCurrentCohortFilters();
  const prevCohortFilters = usePrevious(currentCohortFilters);
  const prevEnumValues = usePrevious(enumValues);

  useEffect(() => {
    if (
      !facet ||
      !isEqual(prevCohortFilters, currentCohortFilters) ||
      !isEqual(prevEnumValues, enumValues)
    ) {
      coreDispatch(
        fetchFacetByNameGQL({
          field: field,
          docType: docType,
          index: indexType,
          caseFilterSelector: selectCurrentCohortFilters,
        }),
      );
    }
  }, [
    coreDispatch,
    facet,
    field,
    docType,
    indexType,
    prevCohortFilters,
    prevEnumValues,
    enumValues,
    currentCohortFilters,
  ]);

  return {
    data: facet?.buckets,
    error: facet?.error,
    isUninitialized: facet === undefined,
    isFetching: facet?.status === "pending",
    isSuccess: facet?.status === "fulfilled",
    isError: facet?.status === "rejected",
  };
};

/**
 * Fetch multiple enum facets via a single API call. Does not return any value but initiates the
 * fetch for each field.
 * @param docType - "cases", "files", etc
 * @param indexType - "explore" | "repository"
 * @param fields - list of fields.
 */
export const useEnumFacets = (
  fields: ReadonlyArray<string>,
  { docType, indexType }: { docType: GQLDocType; indexType: GQLIndexType },
): void => {
  const facet: ReadonlyArray<FacetBuckets> = useCoreSelector((state) =>
    selectMultipleFacetsByDocTypeAndField(state, docType, fields),
  );
  const coreDispatch = useCoreDispatch();

  const enumValues = useEnumFiltersByNames(fields);
  const currentCohortFilters = useCohortFacetFilter();
  const prevCohortFilters = usePrevious(currentCohortFilters);
  const prevEnumValues = usePrevious(enumValues);
  const prevFilterLength = usePrevious(facet.length);

  useEffect(() => {
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    if (
      facet.length > 0 &&
      (prevFilterLength != facet.length ||
        !isEqual(prevCohortFilters, currentCohortFilters) ||
        !isEqual(prevEnumValues, enumValues))
    ) {
      coreDispatch(
        fetchFacetByNameGQL({
          field: fields,
          docType: docType,
          index: indexType,
          caseFilterSelector: selectCurrentCohortFilters,
        }),
      );
    }
  }, [
    coreDispatch,
    facet,
    fields,
    docType,
    indexType,
    prevCohortFilters,
    prevEnumValues,
    enumValues,
    currentCohortFilters,
    prevFilterLength,
  ]);
};

export const useAllEnumFacets = () => {
  const tabsConfig = useCoreSelector((state) =>
    selectCohortBuilderConfig(state),
  );
  const facets = useCoreSelector((state) => selectFacetDefinition(state));

  const coreDispatch = useCoreDispatch();
  const currentCohortFilters = useCurrentCohortFilters();

  useDeepCompareEffect(() => {
    Object.entries(tabsConfig).map(([, tabEntry]) => {
      const facetList = getFacetInfo(tabEntry.facets, {
        ...(facets.data || {}),
      });

      const enumFacets = facetList
        .filter((x) => x.facet_type === "enum")
        .map((x) => x.field);

      if (enumFacets.length > 0) {
        coreDispatch(
          fetchFacetByNameGQL({
            field: enumFacets,
            docType: tabEntry.queryOptions.docType as GQLDocType,
            index: tabEntry.queryOptions.indexType as GQLIndexType,
            caseFilterSelector: selectCurrentCohortFilters,
          }),
        );
      }
    });
  }, [tabsConfig, facets, currentCohortFilters, coreDispatch]);
};

type UpdateEnumFiltersFunc = (
  dispatch: ThunkDispatch<any, undefined, UnknownAction>,
  field: string,
  enumerationFilters: EnumOperandValue,
) => void;

/**
 * Adds an enumeration filter to cohort filters
 * @param dispatch - CoreDispatch instance
 * @param enumerationFilters - values to update
 * @param field - field to update
 */
export const updateEnumFilters: UpdateEnumFiltersFunc = (
  dispatch: ThunkDispatch<any, undefined, UnknownAction>,
  field: string,
  enumerationFilters: EnumOperandValue,
) => {
  // undefined just return
  if (enumerationFilters === undefined) return;
  if (enumerationFilters.length > 0) {
    dispatch(
      updateCohortFilter({
        field: field,
        operation: {
          operator: "includes",
          field: field,
          operands: enumerationFilters,
        },
      }),
    );
  } else {
    // completely remove the field
    dispatch(removeCohortFilter(field));
  }
};

export const useRangeFacet = (
  field: string,
  ranges: ReadonlyArray<NumericFromTo>,
  { docType, indexType }: { docType: GQLDocType; indexType: GQLIndexType },
  overrideCohortFilters?: GqlOperation,
) => {
  const coreDispatch = useCoreDispatch();
  const facet: FacetBuckets = useCoreSelector((state) =>
    selectRangeFacetByField(state, field),
  );

  const cohortFilters = useCohortFacetFilter();
  const prevFilters = usePrevious(cohortFilters);
  const prevRanges = usePrevious(ranges);

  const rangeCohortFilters = useMemo(
    () => overrideCohortFilters ?? buildCohortGqlOperator(cohortFilters),
    [overrideCohortFilters, cohortFilters],
  );
  const prevRangeFilters = usePrevious(rangeCohortFilters);

  useEffect(() => {
    if (
      !facet ||
      !isEqual(prevRangeFilters, rangeCohortFilters) ||
      !isEqual(ranges, prevRanges)
    ) {
      coreDispatch(
        fetchFacetContinuousAggregation({
          field: field,
          ranges: ranges,
          docType: docType,
          indexType: indexType,
          overrideFilters: rangeCohortFilters,
        }),
      );
    }
  }, [
    coreDispatch,
    facet,
    field,
    rangeCohortFilters,
    prevFilters,
    prevRangeFilters,
    ranges,
    prevRanges,
    docType,
    indexType,
  ]);

  return {
    data: facet?.buckets,
    error: facet?.error,
    isUninitialized: facet === undefined,
    isFetching: facet?.status === "pending",
    isSuccess: facet?.status === "fulfilled",
    isError: facet?.status === "rejected",
  };
};

// Global Selector for Facet Values
export const useSelectFieldFilter = (field: string): Operation => {
  // get the current filter for this facet
  return useCoreSelector((state) =>
    selectCurrentCohortFiltersByName(state, field),
  );
};

// Update filter hook
export const useUpdateFacetFilter = () => {
  const dispatch = useCoreDispatch();
  // update the filter for this facet
  return (field: string, operation: Operation) => {
    dispatch(updateCohortFilter({ field: field, operation: operation }));
  };
};

// Core ClearFilters hook
export const useClearFilters = () => {
  const dispatch = useCoreDispatch();
  return (field: string) => {
    dispatch(removeCohortFilter(field));
  };
};

export const useEnumFacetValues = (
  field: string,
  queryOptions: { docType: GQLDocType },
) => {
  // facet data is store in core
  const facet: FacetBuckets = useCoreSelector((state) =>
    selectFacetByDocTypeAndField(state, queryOptions.docType, field),
  );

  const enumValues = useCohortFacetFilterByName(field);
  return {
    data: facet?.buckets,
    enumFilters: (enumValues as EnumOperandValue)?.map((x) => x.toString()),
    error: facet?.error,
    isUninitialized: facet === undefined,
    isFetching: facet?.status === "pending",
    isSuccess: facet?.status === "fulfilled",
    isError: facet?.status === "rejected",
  };
};

export const useTotalCounts = ({
  docType,
}: {
  docType: GQLDocType;
}): number => {
  const caseCountName = FacetDocTypeToCountsIndexMap[docType];
  return useCoreSelector((state) =>
    selectTotalCountsByName(state, caseCountName),
  );
};

export const FacetDocTypeToCountsIndexMap = {
  cases: "caseCounts",
  files: "fileCounts",
  genes: "genesCounts",
  ssms: "mutationCounts",
  projects: "projectCounts",
};

export const FacetDocTypeToLabelsMap = {
  cases: "Cases",
  files: "Files",
  genes: "Genes",
  ssms: "Mutations",
  projects: "Projects",
};

export const useLocalFilters = (
  field: string,
  docType: GQLDocType,
  selectFieldEnumValues: (field: string) => OperandValue,
  selectLocalFilters: () => FilterSet,
  facetSelector: (state: CoreState, field: string) => FacetBuckets,
) => {
  const coreDispatch = useCoreDispatch();

  const facet: FacetBuckets = useCoreSelector((state) =>
    facetSelector(state, field),
  ); // Facet data is always cached in the coreState

  const enumValues = selectFieldEnumValues(field);
  const localFilters = selectLocalFilters();
  const prevLocalFilters = usePrevious(localFilters);
  const prevEnumValues = usePrevious(enumValues);

  useEffect(() => {
    if (
      !facet ||
      !isEqual(prevLocalFilters, localFilters) ||
      !isEqual(prevEnumValues, enumValues)
    ) {
      coreDispatch(
        fetchFacetByNameGQL({
          field: field,
          docType: docType,
          localFilters,
        }),
      );
    }
  }, [
    coreDispatch,
    facet,
    field,
    localFilters,
    docType,
    prevLocalFilters,
    prevEnumValues,
    enumValues,
  ]);

  return {
    data: facet?.buckets,
    enumFilters: (enumValues as EnumOperandValue)?.map((x) => x.toString()),
    error: facet?.error,
    isUninitialized: facet === undefined,
    isFetching: facet?.status === "pending",
    isSuccess: facet?.status === "fulfilled",
    isError: facet?.status === "rejected",
  };
};
