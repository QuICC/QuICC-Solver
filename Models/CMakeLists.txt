#
# Register Models
#
#
# A model can be manually added as a folder under QuICC/Models and added
# to the list QUICC_MODEL
#
# Otherwise the models added to the list QUICC_MODEL are looked up in the model repository.
#
#
#
# The following variables affect the behaviour:
#
# QUICC_MODEL
#    List of models to be built (default: )
#    Note that the list is treated as case insensitive (due to GitHub limitation)
#    The case-ness of the first model registration is retained for the path
#
# QUICC_MODEL_GIT_REPO_BASE
#    Git repository to be used instead of default.
#    default: ${QUICC_GITHUB_PREFIX}QuICC/Model-
#    Note that the full repository must be called `<base repo><model>`
#
# QUICC_MODEL_GIT_BRANCH
#    Git repository branch to be used (default: main)
#
# QUICC_MODEL_GIT_REPO_BASE_<model>
#    Git repository to be used for a specific model (default: QUICC_MODEL_GIT_REPO_BASE for unknown models)
#
# QUICC_MODEL_GIT_BRANCH_<model>
#    Git repository branch to be used for a specific model (default: QUICC_MODEL_GIT_BRANCH)
#

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake.d")

set(QUICC_MODEL_PATH "Models")
# Teststuite
include(setup/ModelTestSuite)

include(FetchRepository)

# Default branch for main models
quicc_set_default_git_branch(BoussinesqSphereModifiedTaylor BRANCH "v1.7.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereInertialessDynamo BRANCH "v1.1.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereRDDDynamo BRANCH "v1.0.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereDynamo BRANCH "v2.1.1" TYPE Model)
quicc_set_default_git_branch(BoussinesqShellDynamo BRANCH "v1.14.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereITO BRANCH "v1.1.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqShellITO BRANCH "v1.3.1" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereRDDC BRANCH "v1.1.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereRTC BRANCH "v2.1.1" TYPE Model)
quicc_set_default_git_branch(BoussinesqShellRTC BRANCH "v1.16.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqSphereTC BRANCH "v2.1.1" TYPE Model)
quicc_set_default_git_branch(BoussinesqShellTC BRANCH "v1.13.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqPlaneRBC BRANCH "v1.3.0" TYPE Model)
quicc_set_default_git_branch(BoussinesqPlaneRRBC BRANCH "v1.2.0" TYPE Model)

set(QUICC_MODEL CACHE STRING "List of Models to build.")
message(DEBUG "QUICC_MODEL: ${QUICC_MODEL}")

foreach(_model ${QUICC_MODEL})
  quicc_fetch_repository(${_model} TYPE Model)
endforeach()
