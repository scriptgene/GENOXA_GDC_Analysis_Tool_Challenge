import { createAsyncThunk, createSlice } from "@reduxjs/toolkit";
import { CoreDispatch } from "../../store";
import { CoreState } from "../../reducers";
import { graphqlAPI, GraphQLApiResponse } from "../gdcapi/gdcgraphql";
import { selectCurrentCohortGqlFilters } from "../cohort";
import { convertFacetNameToGQL } from "./facetApiGQL";
import { FacetBuckets, GQLIndexType, GQLDocType } from "./types";
import { RangeBuckets, processRangeResults } from "./continuousAggregationApi";
import { GqlOperation, NumericFromTo } from "../gdcapi/filters";

export const buildContinuousAggregationRangeOnlyQuery = (
  field: string,
  itemType: GQLDocType,
  indexType: GQLIndexType,
  alias?: string,
): string => {
  const queriedFacet =
    alias !== undefined
      ? `${convertFacetNameToGQL(alias)} : ${convertFacetNameToGQL(field)}`
      : convertFacetNameToGQL(field);

  return `
  query ContinuousAggregationQuery($caseFilters: FiltersArgument, $filters: FiltersArgument, $filters2: FiltersArgument) {
  viewer {
    ${indexType} {
      ${itemType} {
        aggregations(case_filters: $caseFilters, filters: $filters) {
          ${queriedFacet} {
           stats {
                Min : min
                Max: max
                Mean: avg
                SD: std_deviation
                count
            }
            range(ranges: $filters2) {
              buckets {
                doc_count
                key
              }
            }
          }
        }
      }
    }
  }
}
`;
};

export interface FetchContinuousAggregationProps {
  readonly field: string;
  readonly ranges: ReadonlyArray<NumericFromTo>;
  readonly docType?: GQLDocType;
  readonly indexType?: GQLIndexType;
  readonly overrideFilters?: GqlOperation;
}

export const fetchFacetContinuousAggregation = createAsyncThunk<
  GraphQLApiResponse<RangeBuckets>,
  FetchContinuousAggregationProps,
  { dispatch: CoreDispatch; state: CoreState }
>(
  "facet/fetchFacetContinuousAggregation",
  async (
    {
      field,
      ranges,
      docType = "cases" as GQLDocType,
      indexType = "explore" as GQLIndexType,
      overrideFilters = undefined,
    },
    thunkAPI,
  ) => {
    const filters = selectCurrentCohortGqlFilters(thunkAPI.getState());
    const adjField = field.includes(`${docType}.`)
      ? field.replace(`${docType}.`, "")
      : field;
    const queryGQL = buildContinuousAggregationRangeOnlyQuery(
      adjField,
      docType,
      indexType,
      field,
    );
    const filtersGQL = {
      caseFilters: overrideFilters ?? filters ?? {},
      filters: {},
      filters2: { op: "range", content: [{ ranges: ranges }] },
    };

    return await graphqlAPI(queryGQL, filtersGQL);
  },
);

const initialState: Record<string, FacetBuckets> = {};

const rangeFacetAggregation = createSlice({
  name: "facet/rangeFacetAggregation",
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(fetchFacetContinuousAggregation.fulfilled, (state, action) => {
        const response = action.payload;
        const index = action.meta.arg.indexType ?? "explore";
        const itemType = action.meta.arg.docType ?? "cases";
        if (response.errors && Object.keys(response.errors).length > 0) {
          state[action.meta.arg.field] = {
            status: "rejected",
            error: JSON.stringify(response.errors),
          };
        } else {
          const aggregations =
            Object(response).data.viewer[index][itemType].aggregations;
          aggregations &&
            processRangeResults(action.meta.requestId, aggregations, state);
        }
      })
      .addCase(fetchFacetContinuousAggregation.pending, (state, action) => {
        const field = action.meta.arg.field;
        state[field] = {
          status: "pending",
          requestId: action.meta.requestId,
        };
      })
      .addCase(fetchFacetContinuousAggregation.rejected, (state, action) => {
        const field = action.meta.arg.field;
        state[field] = {
          status: "rejected",
        };
      });
  },
});

export const rangeFacetsReducer = rangeFacetAggregation.reducer;

export const selectRangeFacets = (
  state: CoreState,
): Record<string, FacetBuckets> => state.facetsGQL.ranges;

export const selectRangeFacetByField = (
  state: CoreState,
  field: string,
): FacetBuckets => state.facetsGQL.ranges[field];
