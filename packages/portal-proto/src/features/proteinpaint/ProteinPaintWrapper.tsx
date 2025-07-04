import { useRef, useCallback, useState, FC } from "react";
import { useDeepCompareEffect } from "use-deep-compare";
import { bindProteinPaint } from "@sjcrh/proteinpaint-client";
import { useIsDemoApp } from "@/hooks/useIsDemoApp";
import {
  useCoreSelector,
  selectCurrentCohortFilters,
  FilterSet,
  PROTEINPAINT_API,
  useFetchUserDetailsQuery,
  useCoreDispatch,
  buildCohortGqlOperator,
  useCreateCaseSetFromValuesMutation,
} from "@gff/core";
import { isEqual, cloneDeep } from "lodash";
import { DemoText } from "@/components/tailwindComponents";
import { SaveCohortModal } from "@gff/portal-components";
import {
  SelectSamples,
  SelectSamplesCallBackArg,
  SelectSamplesCallback,
} from "./sjpp-types";
import { cohortActionsHooks } from "../cohortBuilder/CohortManager/cohortActionHooks";
import { INVALID_COHORT_NAMES } from "../cohortBuilder/utils";

const basepath = PROTEINPAINT_API;

interface PpProps {
  basepath?: string;
  geneId?: string;
  gene2canonicalisoform?: string;
  ssm_id?: string;
  mds3_ssm2canonicalisoform?: mds3_isoform;
  geneSearch4GDCmds3?: boolean;
  hardcodeCnvOnly?: boolean;
}

export const ProteinPaintWrapper: FC<PpProps> = (props: PpProps) => {
  const isDemoMode = useIsDemoApp();
  const currentCohort = useCoreSelector(selectCurrentCohortFilters);
  const filter0 = isDemoMode ? null : buildCohortGqlOperator(currentCohort);
  const userDetails = useFetchUserDetailsQuery();

  // to track reusable instance for mds3 skewer track
  const prevArg = useRef<any>({});
  const coreDispatch = useCoreDispatch();
  const [showSaveCohort, setShowSaveCohort] = useState(false);
  const [createSet, response] = useCreateCaseSetFromValuesMutation();
  const [newCohortFilters, setNewCohortFilters] =
    useState<FilterSet>(undefined);

  const callback = useCallback<SelectSamplesCallback>(
    (arg: SelectSamplesCallBackArg) => {
      const cases = arg.samples.map((d) => d["cases.case_id"]);
      if (cases.length > 1) {
        createSet({ values: cases, intent: "portal", set_type: "frozen" });
      } else {
        setNewCohortFilters({
          mode: "and",
          root: {
            "cases.case_id": {
              operator: "includes",
              field: "cases.case_id",
              operands: cases,
            },
          },
        });
        setShowSaveCohort(true);
      }
    },
    [createSet],
  );

  // a set for the new cohort is created, now show the save cohort modal
  useDeepCompareEffect(() => {
    if (response.isSuccess) {
      const filters: FilterSet = {
        mode: "and",
        root: {
          "cases.case_id": {
            operator: "includes",
            field: "cases.case_id",
            operands: [`set_id:${response.data}`],
          },
        },
      };
      setNewCohortFilters(filters);
      setShowSaveCohort(true);
    }
  }, [response.isSuccess, coreDispatch, response.data]);

  useDeepCompareEffect(
    () => {
      const rootElem = divRef.current as HTMLElement;
      const data = getLollipopTrack(props, filter0, callback);
      if (!data) return;
      if (isDemoMode) {
        data.geneSymbol = props.hardcodeCnvOnly
          ? "chr8:127682515-127792250"
          : "MYC";
      }
      // compare the argument to runpp to avoid unnecessary render
      if ((data || prevArg.current) && isEqual(prevArg.current, data)) return;
      prevArg.current = data;

      const toolContainer = rootElem.parentNode.parentNode
        .parentNode as HTMLElement;
      toolContainer.style.backgroundColor = "#fff";

      const arg = Object.assign(
        { holder: rootElem, noheader: true, nobox: true },
        cloneDeep(data),
      ) as Mds3Arg;

      // bindProteinPaint() handles rapid update requests/race condition,
      // so no need to include debouncing and promise code in this wrapper
      // TODO: will revert to using runproteinpaint() once these advanced capabilities
      // are merged into it
      bindProteinPaint({
        rootElem,
        initArgs: arg,
        updateArgs: arg,
        isStale() {
          // new data has replaced this one, will prevent unnecessary render
          // in case of race condition
          return prevArg.current != data;
        },
      });
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [
      props.gene2canonicalisoform,
      props.mds3_ssm2canonicalisoform,
      props.geneSearch4GDCmds3,
      isDemoMode,
      filter0,
      userDetails,
    ],
  );

  const divRef = useRef();
  const demoText = props.hardcodeCnvOnly
    ? "Demo showing MYC CNV segments from all GDC cases"
    : "Demo showing MYC mutations from all GDC cases.";
  return (
    <div>
      {isDemoMode && <DemoText>{demoText}</DemoText>}
      <div
        ref={divRef}
        className="sjpp-wrapper-root-div"
        //userDetails={userDetails}
      />

      <SaveCohortModal // Show the modal, create a saved cohort when save button is clicked
        opened={showSaveCohort}
        onClose={() => setShowSaveCohort(false)}
        filters={newCohortFilters}
        hooks={cohortActionsHooks}
        invalidCohortNames={INVALID_COHORT_NAMES}
      />
    </div>
  );
};

interface Mds3Arg {
  holder?: HTMLElement;
  noheader?: boolean;
  nobox?: boolean;
  hide_dsHandles?: boolean;
  host: string;
  gene2canonicalisoform?: string;
  mds3_ssm2canonicalisoform?: mds3_isoform;
  geneSearch4GDCmds3?:
    | boolean
    | {
        hardcodeCnvOnly?: boolean;
      };
  geneSymbol?: string;
  tracks?: Track[];
  filter0?: FilterSet;
  allow2selectSamples?: SelectSamples;
}

interface Track {
  type: string;
  dslabel: string;
  filter0: FilterSet;
  allow2selectSamples?: SelectSamples;
  hardcodeCnvOnly?: boolean;
}

interface mds3_isoform {
  ssm_id: string;
  dslabel: string;
}

function getLollipopTrack(
  props: PpProps,
  filter0: any,
  callback: SelectSamplesCallback,
) {
  const arg: Mds3Arg = {
    // host in gdc is just a relative url path,
    // using the same domain as the GDC portal where PP is embedded
    host: props.basepath || (basepath as string),
    geneSearch4GDCmds3: {
      hardcodeCnvOnly: props.hardcodeCnvOnly,
    },
    filter0,
    allow2selectSamples: {
      buttonText: "Create Cohort",
      attributes: [{ from: "sample_id", to: "cases.case_id", convert: true }],
      callback,
    },
  };

  if (props.hardcodeCnvOnly) {
    arg.geneSearch4GDCmds3 = {
      hardcodeCnvOnly: true,
    };
  } else if (props.geneId) {
    arg.gene2canonicalisoform = props.geneId;
  } else if (props.ssm_id) {
    arg.mds3_ssm2canonicalisoform = {
      dslabel: "GDC",
      ssm_id: props.ssm_id,
    };
  } else {
    arg.geneSearch4GDCmds3 = true;
  }

  return arg;
}
