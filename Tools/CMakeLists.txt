#
# Register Tools
#
#
# A tool can be manually added as a folder under QuICC/Tools and added
# to the list QUICC_TOOL
#
# Otherwise the models added to the list QUICC_TOOL are looked up in the tool repository.
#
#
#
# The following variables affect the behaviour:
#
# QUICC_TOOL
#    List of tools to be built (default: Simulation)
#    Note that the list is treated as case insensitive (due to GitHub limitation)
#    The case-ness of the first tool registration is retained for the path
#
# QUICC_TOOL_GIT_REPO_BASE
#    Git repository to be used instead of default.
#    default: ${QUICC_GITHUB_PREFIX}QuICC/Tool-
#    Note that the full repository must be called `<base repo><tool>`
#
# QUICC_TOOL_GIT_BRANCH
#    Git repository branch to be used (default: main)
#
# QUICC_TOOL_GIT_REPO_BASE_<tool>
#    Git repository to be used for a specific tool (default: QUICC_TOOL_GIT_REPO_BASE for unknown models)
#
# QUICC_TOOL_GIT_BRANCH_<tool>
#    Git repository branch to be used for a specific tool (default: QUICC_TOOL_GIT_BRANCH)
#

include(FetchRepository)

# Default branch for main tools
quicc_set_default_git_branch("Simulation" BRANCH "dev" TYPE "Tool")
quicc_set_default_git_branch("Stability" BRANCH "dev" TYPE "Tool")

set(QUICC_TOOL "Simulation" CACHE STRING "List of Tools to build.")
message(DEBUG "QUICC_TOOL: ${QUICC_TOOL}")

foreach(_tool ${QUICC_TOOL})
  quicc_fetch_repository(${_tool} TYPE Tool MAKE_AVAILABLE ON)
endforeach()
