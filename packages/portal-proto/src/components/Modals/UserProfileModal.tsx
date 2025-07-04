import { CheckIcon } from "@/utils/icons";
import { useFetchUserDetailsQuery } from "@gff/core";
import { Text } from "@mantine/core";
import { ScrollableTableWithFixedHeader } from "../ScrollableTableWithFixedHeader/ScrollableTableWithFixedHeader";
import { BaseModal } from "./BaseModal";
import { SessionExpireModal } from "./SessionExpireModal";

export const UserProfileModal = ({
  openModal,
}: {
  openModal: boolean;
}): JSX.Element => {
  const { data: userInfo } = useFetchUserDetailsQuery();

  if (userInfo?.status === 401) {
    return <SessionExpireModal openModal={openModal} />;
  }

  const {
    projects: { gdc_ids },
    username,
  } = userInfo?.data || { username: undefined, projects: {} };

  // get the unique permission properties
  const allPermissionValues = Array.from(
    new Set(
      gdc_ids &&
        Object.keys(gdc_ids).reduce(
          (acc, projectId) => [...acc, ...gdc_ids[projectId]],
          [],
        ),
    ),
  );

  const data =
    gdc_ids &&
    Object.keys(gdc_ids).map((projectId) => ({
      projectId,
      ...allPermissionValues.reduce(
        (acc, v) => ({
          ...acc,
          [v]: gdc_ids[projectId].includes(v) ? (
            <span>
              <CheckIcon
                aria-label="permission given"
                aria-hidden
                key={projectId}
              />
            </span>
          ) : (
            <span aria-label="no permission"></span>
          ),
        }),
        [],
      ),
    }));

  const headings = [
    "Project ID",
    ...allPermissionValues.map((v) => {
      if (v === "_member_") {
        return "member";
      }
      return v;
    }),
  ];

  return (
    <BaseModal
      title={
        <Text size="lg" className="font-medium">{`Username: ${username}`}</Text>
      }
      openModal={openModal}
      size="60%"
      buttons={[{ title: "Close", dataTestId: "button-close" }]}
    >
      <div className={`${!data ? "py-4" : "py-2"}`}>
        {data?.length > 0 ? (
          <ScrollableTableWithFixedHeader
            customDataTestID="table-user-profile-access"
            tableData={{
              headers: headings,
              tableRows: data,
            }}
            maxRowsBeforeScroll={10}
          />
        ) : (
          <div data-testid="warningText" className="mt-4">
            <Text className="mb-4 font-content">
              You do not have any access to controlled access data for projects
              available in the GDC Data Portal.
            </Text>
            <Text className="font-content">
              For instructions on{" "}
              <a
                href="https://gdc.cancer.gov/access-data/obtaining-access-controlled-data"
                target="_blank"
                rel="noreferrer"
                className="text-utility-link underline hover:text-accent-cool-darkest"
              >
                how to apply for access to controlled data
              </a>
              , please visit our documentation on how to apply for access
              through dbGAP.
            </Text>
          </div>
        )}
      </div>
    </BaseModal>
  );
};
