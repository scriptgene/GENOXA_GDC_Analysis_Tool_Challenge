import { useRef, FC } from "react";
import { useDeepCompareEffect } from "use-deep-compare";
import { runproteinpaint } from "@sjcrh/proteinpaint-client";
import {
  useCoreSelector,
  selectCurrentCohortFilters,
  PROTEINPAINT_API,
  useFetchUserDetailsQuery,
  buildCohortGqlOperator,
} from "@gff/core";
import { isEqual, cloneDeep } from "lodash";

const basepath = PROTEINPAINT_API;

interface PpProps {
  basepath?: string;
}

export const CohortLevelMafWrapper: FC<PpProps> = (props: PpProps) => {
  // to track reusable instance for mds3 skewer track
  const ppRef = useRef<PpApi>();
  const prevArg = useRef<any>({});
  const currentCohort = useCoreSelector(selectCurrentCohortFilters);
  const filter0 = buildCohortGqlOperator(currentCohort);
  const userDetails = useFetchUserDetailsQuery();

  useDeepCompareEffect(
    () => {
      const rootElem = divRef.current as HTMLElement;
      const data = getMafUi(props, filter0);
      if (!data) return;
      // compare the argument to runpp to avoid unnecessary render
      if ((data || prevArg.current) && isEqual(prevArg.current, data)) return;
      prevArg.current = data;

      const toolContainer = rootElem.parentNode.parentNode
        .parentNode as HTMLElement;
      toolContainer.style.backgroundColor = "#fff";

      const arg = Object.assign(
        { holder: rootElem, noheader: true, nobox: true, hide_dsHandles: true },
        cloneDeep(data),
      ) as MafUiArg;

      if (ppRef.current) {
        ppRef.current.update(arg);
      } else {
        const pp_holder = rootElem.querySelector(".sja_root_holder");
        if (pp_holder) pp_holder.remove();
        runproteinpaint(arg).then((pp) => {
          ppRef.current = pp;
        });
      }
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [filter0, userDetails],
  );

  const divRef = useRef();
  return (
    <div>
      <div
        ref={divRef}
        className="sjpp-wrapper-root-div"
        //userDetails={userDetails}
      />
    </div>
  );
};

interface MafUiArg {
  holder?: HTMLElement;
  host: string;
  noheader?: boolean;
  nobox?: boolean;
  hide_dsHandles?: boolean;
  launchGdcMaf: true;
  genome?: string;
  filter0: any;
}

interface PpApi {
  update(arg: any): null;
}

function getMafUi(props: PpProps, filter0: any) {
  const arg: MafUiArg = {
    // host in gdc is just a relative url path,
    // using the same domain as the GDC portal where PP is embedded
    host: props.basepath || (basepath as string),
    filter0: filter0 || null,
    launchGdcMaf: true,
  };
  return arg;
}
