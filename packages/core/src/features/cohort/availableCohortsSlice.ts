import {
  createEntityAdapter,
  createSlice,
  EntityState,
  nanoid,
  type PayloadAction,
  createAsyncThunk,
  ThunkAction,
  UnknownAction,
} from "@reduxjs/toolkit";
import { CoreState } from "../../reducers";
import { buildCohortGqlOperator, FilterSet } from "./filters";
import { GqlOperation, Operation } from "../gdcapi/filters";
import { CoreDataSelectorResponse, DataStatus } from "../../dataAccess";
import { CoreDispatch } from "../../store";
import { useCoreSelector } from "../../hooks";
import { defaultCohortNameGenerator } from "./utils";
import { CountsData, CountsDataAndStatus, NullCountsData } from "./type";
import { graphqlAPI, GraphQLApiResponse } from "../gdcapi/gdcgraphql";
import { UnknownJson } from "../gdcapi/types";

export const UNSAVED_COHORT_NAME = "Unsaved_Cohort";

export interface CaseSetDataAndStatus {
  readonly status: DataStatus; // status of create caseSet
  readonly error?: string; // any error message
  readonly caseSetIds?: Record<string, string>; // prefix mapped caseSetIds
}

/**
 * A Cohort is a collection of filters that can be used to query the GDC API.
 * The cohort interface is used to manage the cohort state in the redux-toolkit entity adapter.
 * @see https://redux-toolkit.js.org/api/createEntityAdapter
 *
 * @property id - the id of the cohort
 * @property name - the name of the cohort
 * @property filters - the filters for the cohort
 * @property caseSet - the caseSet for the cohort
 * @property sets - the sets for the cohort
 * @property modified - flag indicating if the cohort has been modified
 * @property modified_datetime - the last time the cohort was modified
 * @property saved - flag indicating if the cohort has been saved
 * @category Cohort
 */

export interface Cohort {
  readonly id: string;
  readonly name: string;
  readonly filters: FilterSet; // active filters for cohort
  readonly caseSet: CaseSetDataAndStatus; // case ids for frozen cohorts
  readonly modified?: boolean; // flag which is set to true if modified and unsaved
  readonly modified_datetime: string; // last time cohort was modified
  readonly saved?: boolean; // flag indicating if cohort has been saved.
  readonly counts: CountsDataAndStatus; //case, file, etc. counts of a cohort
  readonly unsavedCohortId?: string; // the id before cohort is saved, used for apps where the cohort can be saved while the user is operating the app
}

const cohortsAdapter = createEntityAdapter<Cohort, CohortId>({
  sortComparer: (a, b) => {
    if (a.modified_datetime <= b.modified_datetime) return 1;
    else return -1;
  },
  selectId: (cohort: Cohort): CohortId => cohort.id,
});

export interface CurrentCohortState {
  currentCohort: string | undefined;
  isLoggedIn?: boolean | undefined;
}

const emptyInitialState = cohortsAdapter.getInitialState<CurrentCohortState>({
  currentCohort: undefined,
  isLoggedIn: undefined, // isLoggedIn is used to trigger an api call to fetch fresh data when user logs in or out of the app
});

interface UpdateFilterParams {
  field: string;
  operation: Operation;
}

export const createCohortName = (postfix: string): string => {
  return `Custom Cohort ${postfix}`;
};

export type CohortId = string;

export const createCohortId = (): CohortId => nanoid();

const newCohort = ({
  filters = { mode: "and", root: {} },
  customName,
}: {
  filters?: FilterSet;
  customName?: string;
}): Cohort => {
  const ts = new Date();
  const newName = customName ?? defaultCohortNameGenerator();

  const newId = createCohortId();
  return {
    name: newName,
    id: newId,
    caseSet: {
      caseSetIds: undefined,
      status: "uninitialized" as DataStatus,
    },
    filters: filters,
    modified: false,
    saved: false,
    modified_datetime: ts.toISOString(),
    counts: {
      ...NullCountsData,
    },
  };
};

const getCurrentCohort = (
  state: EntityState<Cohort, CohortId> & CurrentCohortState,
): CohortId => {
  if (state.currentCohort) {
    return state.currentCohort;
  }

  const unsavedCohort = newCohort({ customName: UNSAVED_COHORT_NAME });
  return unsavedCohort.id;
};

const getCurrentCohortFromCoreState = (state: CoreState): CohortId => {
  if (state.cohort.availableCohorts.currentCohort) {
    return state.cohort.availableCohorts.currentCohort;
  }
  const unsavedCohort = newCohort({ customName: UNSAVED_COHORT_NAME });
  return unsavedCohort.id;
};

const checkForUnsavedCohorts = (
  state: EntityState<Cohort, CohortId> & CurrentCohortState,
  replace: boolean,
) => {
  const selector = cohortsAdapter.getSelectors();
  const unsavedCohorts = selector
    .selectAll(state)
    .filter((cohort) => !cohort.saved);
  if (unsavedCohorts.length > 0) {
    if (replace) {
      cohortsAdapter.removeMany(
        state,
        unsavedCohorts.map((c) => c.id),
      );
    } else {
      throw "There is a limit of one unsaved cohort at a time for a user. Please create a saved cohort or replace the current unsaved cohort";
    }
  }
};

interface NewUnsavedCohortParams {
  filters: FilterSet; // set the filters for the new cohort
  name: string; // set the name for the new cohort
  replace?: boolean; // Replace the current unsaved cohort if there is one
}

interface CopyCohortParams {
  sourceId: string;
  destId: string;
}

/**
 * A Cohort Management Slice which allow cohort to be created and updated.
 * this uses redux-toolkit entity adapter to manage the cohorts
 * Because it is an entity adapter, the state contains an array of id (string)
 * and a Dictionary of Cohort objects. There are two additional members:
 *  - currentCohortId: which is used to identify the "current" or active cohort
 *
 * The slice exports the following actions:
 * - setCohortList() - set saved cohort to the adapter that comes from the server
 * - addNewDefaultUnsavedCohort - create a an instance of the default unsaved cohort
 * - addNewSavedCohort - add a saved cohort
 * - addNewUnsavedCohort - create a new unsaved cohort with the passed filters
 * - copyToSavedCohort - create a copy of the cohort with sourceId to a new cohort with destId
 * - updateCohortName(name:string): changes the current cohort's name
 * - updateCohortFilter(filters: FilterSet): update the filters for this cohort
 * - removeCohortFilter(filter:string): removes the filter from the cohort
 * - clearCohortFilters(): removes all the filters by setting them to the default all GDC state
 * - setCurrentCohortId(id:string): set the id of the current cohort, used to switch between cohorts
 * - clearCaseSet(): resets the caseSet member to all GDC
 * - removeCohort(): removes the current cohort
 * - addNewCohortGroups(): adds groups of filters to the current cohort
 * - removeCohortGroup(): removes a group of filters from the current cohort
 * @category Cohort
 */
const slice = createSlice({
  name: "cohort/availableCohorts",
  initialState: emptyInitialState,
  reducers: {
    /** @hidden */
    setCohortList: (state, action: PayloadAction<Cohort[]>) => {
      // TODO: Behavior TBD - https://jira.opensciencedatacloud.org/browse/PEAR-762
      // When the user deletes context id from their cookies
      // All the cohorts that was previously were saved or unsaved should be removed from the adapter
      if (!action.payload) {
        cohortsAdapter.removeMany(state, state.ids);
      } else {
        cohortsAdapter.upsertMany(state, [...action.payload] as Cohort[]);
      }
    },
    addNewDefaultUnsavedCohort: (state) => {
      const cohort = newCohort({
        customName: UNSAVED_COHORT_NAME,
      });
      checkForUnsavedCohorts(state, false);
      cohortsAdapter.addOne(state, cohort);
      state.currentCohort = cohort.id;
    },
    addNewSavedCohort: (state, action: PayloadAction<Cohort>) => {
      cohortsAdapter.setOne(state, {
        ...action.payload,
        modified: false,
        saved: true,
      });
    },
    addNewUnsavedCohort: (
      state,
      action: PayloadAction<NewUnsavedCohortParams>,
    ) => {
      const cohort = newCohort({
        filters: action.payload.filters,
        customName: action.payload.name,
      });
      checkForUnsavedCohorts(state, action.payload?.replace || false);
      cohortsAdapter.addOne(state, cohort);
      state.currentCohort = cohort.id;
    },
    copyToSavedCohort: (state, action: PayloadAction<CopyCohortParams>) => {
      const sourceCohort = state.entities[action.payload.sourceId];
      if (sourceCohort) {
        const destCohort = {
          ...sourceCohort,
          id: action.payload.destId,
          modified: false,
          saved: true,
          unsavedCohortId: action.payload.sourceId,
          counts: {
            // will need to re-request counts
            ...NullCountsData,
          },
        };
        cohortsAdapter.addOne(state, destCohort);
      }
    },
    updateCohortName: (state, action: PayloadAction<string>) => {
      cohortsAdapter.updateOne(state, {
        id: getCurrentCohort(state),
        changes: { name: action.payload },
      });
    },
    removeCohort: (
      state,
      action: PayloadAction<{
        id?: string;
      }>,
    ) => {
      cohortsAdapter.removeOne(
        state,
        action?.payload?.id || getCurrentCohort(state),
      );

      const selector = cohortsAdapter.getSelectors();
      if (
        selector.selectAll(state).length > 0 &&
        action?.payload.id === undefined
      ) {
        state.currentCohort = selector.selectAll(state)[0].id;
      }
    },
    updateCohortFilter: (state, action: PayloadAction<UpdateFilterParams>) => {
      const filters = {
        mode: "and",
        root: {
          ...state.entities[getCurrentCohort(state)]?.filters.root,
          [action.payload.field]: action.payload.operation,
        },
      };

      cohortsAdapter.updateOne(state, {
        id: getCurrentCohort(state),
        changes: {
          filters: filters,
          modified: true,
          modified_datetime: new Date().toISOString(),
        },
      });
    },
    removeCohortFilter: (state, action: PayloadAction<string>) => {
      const root = state.entities[getCurrentCohort(state)]?.filters.root;
      if (!root) {
        return;
      }
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const { [action.payload]: _a, ...updated } = root;

      cohortsAdapter.updateOne(state, {
        id: getCurrentCohort(state),
        changes: {
          filters: { mode: "and", root: updated },
          modified: true,
          modified_datetime: new Date().toISOString(),
          caseSet: {
            caseSetIds: undefined,
            status: "uninitialized",
          },
        },
      });
    },
    clearCohortFilters: (state) => {
      cohortsAdapter.updateOne(state, {
        id: getCurrentCohort(state),
        changes: {
          filters: { mode: "and", root: {} },
          modified: true,
          modified_datetime: new Date().toISOString(),
          caseSet: {
            caseSetIds: undefined,
            status: "uninitialized",
          },
        },
      });
    },
    discardCohortChanges: (
      state,
      action: PayloadAction<{
        filters: FilterSet | undefined;
        id?: string;
      }>,
    ) => {
      cohortsAdapter.updateOne(state, {
        id: action.payload.id ?? getCurrentCohort(state),
        changes: {
          filters: action.payload.filters || { mode: "and", root: {} },
          modified: false,
          modified_datetime: new Date().toISOString(),
        },
      });
    },
    setCurrentCohortId: (state, action: PayloadAction<string>) => {
      state.currentCohort = action.payload;
    },
    setIsLoggedIn: (state, action: PayloadAction<boolean>) => {
      state.isLoggedIn = action.payload;
    },
    clearCaseSet: (state) => {
      cohortsAdapter.updateOne(state, {
        id: getCurrentCohort(state),
        changes: {
          caseSet: {
            status: "uninitialized",
            caseSetIds: undefined,
          },
        },
      });
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(fetchCohortCaseCounts.fulfilled, (state, action) => {
        if (
          action.meta.requestId !==
          state.entities[action.meta.arg]?.counts.requestId
        ) {
          // ignore if the request id is not the same as the current cohort
          return;
        }

        const response = action.payload;
        if (response.errors && Object.keys(response.errors).length > 0) {
          cohortsAdapter.updateOne(state, {
            id: action.meta.arg,
            changes: {
              counts: { ...NullCountsData, status: "rejected" },
            },
          });
          return;
        }

        // copy the counts for explore and repository
        cohortsAdapter.updateOne(state, {
          id: action.meta.arg,
          changes: {
            counts: {
              caseCount: response.data.viewer.explore.cases.hits.total,
              genesCount: response.data.viewer.explore.genes.hits.total,
              mutationCount: response.data.viewer.explore.ssms.hits.total,
              fileCount: response.data.viewer.repository.files.hits.total,
              status: "fulfilled",
            },
          },
        });
      })
      .addCase(fetchCohortCaseCounts.pending, (state, action) => {
        cohortsAdapter.updateOne(state, {
          id: action.meta.arg,
          changes: {
            counts: {
              caseCount: -1,
              fileCount: -1,
              genesCount: -1,
              mutationCount: -1,
              status: "pending",
              requestId: action.meta.requestId,
            },
          },
        });
      })
      .addCase(fetchCohortCaseCounts.rejected, (state, action) => {
        if (
          action.meta.requestId !==
          state.entities[action.meta.arg]?.counts.requestId
        ) {
          return;
        }
        cohortsAdapter.updateOne(state, {
          id: action.meta.arg,
          changes: {
            counts: {
              caseCount: -1,
              fileCount: -1,
              genesCount: -1,
              mutationCount: -1,
              status: "rejected",
            },
          },
        });
      });
  },
});

export const availableCohortsReducer = slice.reducer;

export const {
  addNewDefaultUnsavedCohort,
  addNewUnsavedCohort,
  removeCohort,
  updateCohortName,
  addNewSavedCohort,
  setCurrentCohortId,
  updateCohortFilter,
  removeCohortFilter,
  clearCohortFilters,
  clearCaseSet,
  setCohortList,
  copyToSavedCohort,
  discardCohortChanges,
  setIsLoggedIn,
} = slice.actions;

/**
 * --------------------------------------------------------------------------
 *  Cohort Selectors
 *  -------------------------------------------------------------------------
 **/

/**
 * Returns the selectors for the cohorts EntityAdapter
 * @param state - the CoreState
 *
 * @hidden
 */
export const cohortSelectors = cohortsAdapter.getSelectors(
  (state: CoreState) => state.cohort.availableCohorts,
);

/**
 * Returns all the cohorts in the state
 * @param state - the CoreState
 *
 * @category Cohort
 * @category Selectors
 */

export const selectAvailableCohorts = (state: CoreState): Cohort[] =>
  cohortSelectors.selectAll(state);

/**
 * Returns the current cohort id
 * @param state - the CoreState
 *
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortId = (state: CoreState): string | undefined =>
  state.cohort?.availableCohorts?.currentCohort;

/**
 * Returns the current cohort login status
 * The case composition of a cohort can be different based on the
 * login status of the user. The login status is attached to the cohort
 * so that the portal and applications are updated are when the login state changes.
 * @param state - the CoreState
 * @hidden
 */
export const selectCohortIsLoggedIn = (state: CoreState): boolean | undefined =>
  state.cohort.availableCohorts.isLoggedIn;

/**
 * Returns if the current cohort is modified
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 * @hidden
 */
export const selectCurrentCohortModified = (
  state: CoreState,
): boolean | undefined => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  return cohort?.modified;
};

/**
 * Returns if the current cohort has been saved
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 * @hidden
 */
export const selectCurrentCohortSaved = (
  state: CoreState,
): boolean | undefined => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  return cohort?.saved;
};

/**
 * Returns the current cohort or undefined if cohort is not found
 * @param state - the CoreState
 * @returns the current cohort or undefined
 *
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohort = (state: CoreState): Cohort | undefined =>
  cohortSelectors.selectById(state, getCurrentCohortFromCoreState(state));

/**
 *  Returns the current cohort name
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortName = (
  state: CoreState,
): string | undefined => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  return cohort?.name;
};

/**
 * Finds a cohort by name and returns it or undefined if not found
 * @param state - the CoreState
 * @param name - the name of the cohort
 * @category Cohort
 * @category Selectors
 */

export const selectAvailableCohortByName = (
  state: CoreState,
  name: string,
): Cohort | undefined =>
  cohortSelectors
    .selectAll(state)
    .find((cohort: Cohort) => cohort.name === name);

/**
 * Returns the cohort's name given the id
 * @param state - the CoreState
 * @param cohortId - the cohort id
 * @category Cohort
 * @category Selectors
 */
export const selectCohortNameById = (
  state: CoreState,
  cohortId: string,
): string | undefined => {
  const cohort = cohortSelectors.selectById(state, cohortId);
  return cohort?.name;
};

/**
 * Returns the cohort's filters given an id
 * This includes the `isLoggedIn` state to ensure fresh data is fetched
 * when the user logs in or out to avoid the use of cached data.
 * @param state - the CoreState
 * @param cohortId - the cohort id
 * @category Cohort
 * @category Selectors
 *
 */
export const selectCohortFilterSetById = (
  state: CoreState,
  cohortId: string,
): FilterSet | undefined => {
  const isLoggedIn = state.cohort.availableCohorts.isLoggedIn;
  const cohort = cohortSelectors.selectById(state, cohortId);
  if (cohort === undefined)
    return { mode: "and", root: {}, isLoggedIn: isLoggedIn };
  return { ...cohort.filters, isLoggedIn: isLoggedIn };
};

/**
 * Returns the currentCohortFilters as a GqlOperation
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortGqlFilters = (
  state: CoreState,
): GqlOperation | undefined => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  return buildCohortGqlOperator(cohort?.filters);
};

/**
 * Public selector of the current Cohort Filters.
 * Returns the current cohort filters as a FilterSet
 * This includes the `isLoggedIn` state to ensure fresh data is fetched
 * when the user logs in or out to avoid the use of cached data.
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortFilters = (state: CoreState): FilterSet => {
  const isLoggedIn = state.cohort.availableCohorts.isLoggedIn;
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  if (cohort === undefined)
    return { mode: "and", root: {}, isLoggedIn: isLoggedIn };
  return { ...cohort.filters, isLoggedIn: isLoggedIn };
};

/**
 * Select a filter by its name from the current cohort. If the filter is not found
 * returns undefined.
 * @param state - the CoreState
 * @param name - the name of the filter
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortFiltersByName = (
  state: CoreState,
  name: string,
): Operation | undefined => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  return cohort?.filters?.root[name];
};

/**
 * Returns the current cohort case count
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortCaseCount = (
  state: CoreState,
): number | undefined =>
  cohortSelectors.selectById(state, getCurrentCohortFromCoreState(state))
    ?.counts.caseCount;

export const selectCurrentCohortFiltersByNames = (
  state: CoreState,
  names: ReadonlyArray<string>,
): Record<string, Operation> => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );

  // TODO replace with memoized selector
  return names.reduce((obj, name) => {
    if (cohort?.filters?.root[name]) obj[name] = cohort?.filters?.root[name];
    return obj;
  }, {} as Record<string, Operation>);
};

/**
 * Returns the current caseSetId filter representing the cohort
 * if the cohort is undefined it returns an empty caseSetIdFilter
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectCurrentCohortCaseSet = (
  state: CoreState,
): CoreDataSelectorResponse<FilterSet> => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  // TODO replace with memoized selector
  if (cohort === undefined || cohort?.caseSet === undefined)
    return {
      data: { mode: "and", root: {} },
      status: "uninitialized",
    };
  return { ...cohort.caseSet };
};

/**
 * Given a cohortId returns the Cohort if found, otherwise returns undefined
 * @param state - the CoreState
 * @param cohortId - the cohort id to return
 * @category Cohort
 * @category Selectors
 */
export const selectCohortByIdOrName = (
  state: CoreState,
  cohortId: string,
  cohortName?: string,
): Cohort | undefined => {
  const id = cohortSelectors.selectById(state, cohortId);

  if (id === undefined) {
    const idByUnsaved = cohortSelectors
      .selectAll(state)
      .find((c) => c?.unsavedCohortId === cohortId);
    if (idByUnsaved) {
      return idByUnsaved;
    }

    return cohortSelectors
      .selectAll(state)
      .find((cohort) => cohort.name == cohortName);
  }

  return id;
};

export const selectMultipleCohortsByIdOrName = (
  state: CoreState,
  cohorts: { id: string; name?: string }[],
) =>
  cohorts
    .map((cohort) => selectCohortByIdOrName(state, cohort.id, cohort.name))
    .filter((cohort) => cohort !== undefined);

/**
 * Returns an array of all the cohorts
 * @param state - the CoreState
 * @category Cohort
 * @category Selectors
 */
export const selectAllCohorts = (state: CoreState): Record<CohortId, Cohort> =>
  cohortSelectors.selectEntities(state);

export const selectCohortCountsResults = (
  state: CoreState,
  cohortId: CohortId = getCurrentCohortFromCoreState(state),
): CoreDataSelectorResponse<CountsData> => {
  return {
    data: cohortSelectors.selectById(state, cohortId)?.counts ?? NullCountsData,
    status:
      cohortSelectors.selectById(state, cohortId)?.counts?.status ??
      NullCountsData.status,
  };
};

export const selectCohortCounts = (
  state: CoreState,
  cohortId: CohortId = getCurrentCohortFromCoreState(state),
): CountsDataAndStatus =>
  cohortSelectors.selectById(state, cohortId)?.counts ?? NullCountsData;

export const selectCohortCountsByName = (
  state: CoreState,
  name: keyof CountsData,
  cohortId: CohortId = getCurrentCohortFromCoreState(state),
): number =>
  cohortSelectors.selectById(state, cohortId)?.counts[name] ??
  NullCountsData[name];

export const selectHasUnsavedCohorts = (state: CoreState): boolean =>
  cohortSelectors.selectAll(state).filter((c) => !c.saved).length >= 1;

export const selectUnsavedCohortName = (state: CoreState): string | undefined =>
  cohortSelectors.selectAll(state).find((c) => !c.saved)?.name;

/**
 * --------------------------------------------------------------------------
 *  Cohort Hooks
 *  -------------------------------------------------------------------------
 **/

/**
 * A hook to get the current cohort filter as a FilterSet
 * @category Cohort
 * @category Hooks
 */
export const useCurrentCohortFilters = (): FilterSet | undefined => {
  return useCoreSelector((state: CoreState) =>
    selectCurrentCohortFilters(state),
  );
};

/**
 * A hook to get the counts of cases and files for the current cohort
 * @category Cohort
 * @category Hooks
 */
export const useCurrentCohortCounts =
  (): CoreDataSelectorResponse<CountsData> => {
    return useCoreSelector((state: CoreState) =>
      selectCohortCountsResults(state),
    );
  };

/**
 * --------------------------------------------------------------------------
 *  Helper Functions
 *  -------------------------------------------------------------------------
 **/

/**
 * A thunk to create a case set when adding filter that require them
 * This primary used to handle gene and ssms applications
 * and is also called from the query expression to handle removing
 * genes and ssms from the expression
 * @param field - the field that requires a case set
 * @param operation - the new filter operation
 * @category Cohort
 * @category Thunks
 * @hidden
 */
export const updateActiveCohortFilter =
  ({
    field,
    operation,
  }: UpdateFilterParams): ThunkAction<
    void,
    CoreState,
    undefined,
    UnknownAction
  > =>
  async (dispatch: CoreDispatch /* getState */) => {
    dispatch(updateCohortFilter({ field, operation }));
  };

/**
 * a thunk to optionally create a caseSet when switching cohorts.
 * Note the assumption if the caseset member has ids then the caseset has previously been created.
 */
export const setActiveCohort =
  (cohortId: string): ThunkAction<void, CoreState, undefined, UnknownAction> =>
  async (dispatch: CoreDispatch /* getState */) => {
    dispatch(setCurrentCohortId(cohortId));
  };

export const discardActiveCohortChanges =
  (
    filters: FilterSet,
  ): ThunkAction<void, CoreState, undefined, UnknownAction> =>
  async (dispatch: CoreDispatch /* getState */) => {
    dispatch(discardCohortChanges({ filters }));
  };

export const setActiveCohortList =
  (cohorts: Cohort[]): ThunkAction<void, CoreState, undefined, UnknownAction> =>
  async (dispatch: CoreDispatch) => {
    // set the list of all cohorts
    dispatch(setCohortList(cohorts));

    // have to request counts for all cohorts loaded from the backend
    cohorts.forEach((cohort) => {
      dispatch(fetchCohortCaseCounts(cohort.id));
    });
  };

export const getCohortFilterForAPI = (cohort: Cohort): FilterSet =>
  cohort.filters;

export const divideFilterSetByPrefix = (
  filters: FilterSet,
  prefixes: string[],
): SplitFilterSet => {
  const results = Object.entries(filters.root).reduce(
    (
      newObj: {
        withPrefix: Record<string, Operation>;
        withoutPrefix: Record<string, Operation>;
      },
      [key, val],
    ) => {
      if (prefixes.some((prefix) => key.startsWith(prefix)))
        newObj.withPrefix[key] = val;
      else newObj.withoutPrefix[key] = val;
      return newObj;
    },
    { withPrefix: {}, withoutPrefix: {} },
  );

  return {
    withPrefix: { mode: "and", root: results.withPrefix },
    withoutPrefix: { mode: "and", root: results.withoutPrefix },
  };
};

interface SplitFilterSet {
  withPrefix: FilterSet;
  withoutPrefix: FilterSet;
}

/**
 * Divides the current cohort into prefix and not prefix. This is
 * used to create caseSets for files and cases
 * @param prefixes - and array of filter prefix to separate on (typically ["genes."])
 * @param state - the CoreState
 */
export const divideCurrentCohortFilterSetFilterByPrefix = (
  state: CoreState,
  prefixes: string[],
): SplitFilterSet => {
  const cohort = cohortSelectors.selectById(
    state,
    getCurrentCohortFromCoreState(state),
  );
  if (cohort === undefined)
    return {
      withPrefix: { mode: "and", root: {} },
      withoutPrefix: { mode: "and", root: {} },
    };

  return divideFilterSetByPrefix(cohort?.filters, prefixes);
};

const CountsGraphQLQuery = `
  query countsQuery($filters: FiltersArgument) {
  viewer {
    repository {
      files {
        hits(case_filters: $filters, first: 0) {
          total
        }
      }
    },
    explore {
      cases {
        hits(case_filters: $filters, first: 0) {
          total
        }
      },
      genes {
        hits(case_filters: $filters, first: 0) {
          total
        }
      },
      ssms {
        hits(case_filters: $filters, first: 0) {
          total
        }
      }
    }
  }
}`;

interface CohortCountsResponse extends GraphQLApiResponse {
  cohortFilters?: FilterSet;
}

export const fetchCohortCaseCounts = createAsyncThunk<
  CohortCountsResponse,
  string,
  { dispatch: CoreDispatch; state: CoreState }
>(
  "cohort/CohortCounts",
  async (cohortId, thunkAPI): Promise<CohortCountsResponse> => {
    const cohortFilters = selectCohortFilterSetById(
      thunkAPI.getState(),
      cohortId,
    );

    const cohortFiltersGQL = buildCohortGqlOperator(cohortFilters);
    const graphQlFilters = {
      filters: cohortFiltersGQL ?? {}, // the cohort filters
    };
    // get the data from the graphql endpoint
    const data = await graphqlAPI<UnknownJson>(
      CountsGraphQLQuery,
      graphQlFilters,
    );
    return {
      ...data,
      cohortFilters, // add the cohort filters to the response to ensure we did not receive stale data
    };
  },
);
