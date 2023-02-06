#pragma once
#include <string>
/**
 * \file
 *  File that provides definitions of all group structures used in algorithms
 */
const std::string GEOM_GROUP_NAME = "/GEOM_GROUP";
const std::string DOSE_GROUP_NAME = "/DOSE_GROUP";
const std::string DAMAGE_GROUP_NAME = "/DAMAGE_GROUP";
const std::string FLOW_GROUP_NAME = "/FLOW_GROUP";
const std::string NUM_LEVELS_DATASET = "/numLevels";
const std::string BOUNDING_BOX_DATASET = "/boundingBox";
const std::string GEOM_DATASET = "/GEOM_ARRAY";
const std::string NODE_DATASET = "/NODE_ARRAY";
const std::string HEALTHY_FLOW_DATASET = "/HEALTHY_FLOW_ARRAY";
const std::string DAMAGED_FLOW_DATASET = "/DAMAGED_FLOW_ARRAY";
const std::string HEALTHY_CONDUCTANCE_DATASET = "/HEALTHY_CONDUCTANCE_ARRAY";
const std::string CONNECTED_VESSELS_DATASET = "/CONNECTED_VESSELS_ARRAY";
const std::string NUM_CONNECTED_VESSELS_DATASET = "/NUM_CONNECTED_VESSELS_ARRAY";
const std::string BOUNDARY_CONDITION_NODE_DATASET = "/BOUNDARY_CONDITION_NODES";
const std::string BOUNDARY_CONDITION_VALUE_DATASET = "/BOUNDARY_CONDITION_VALUE";
const std::string BOUNDARY_CONDITION_TYPE_DATASET = "/BOUNDARY_CONDITION_TYPE";
const std::string NUM_VESSELS_PER_PARTITION_DATASET = "/VESSELS_PER_PARTITION_ARRAY";
const std::string ORGINAL_INDEX_DATASET = "/ORIGINAL_INDEX_ARRAY";
const std::string EDGE_NODE_INDEX_DATASET = "/EDGE_NODE_INDEX_ARRAY";
const std::string EDGE_NODE_HOST_PARTITION_DATASET =
    "/EDGE_NODE_HOST_PARTITION_ARRAY";
const std::string EDGE_NODE_NUM_ASSOCIATED_PARTITIONS_DATASET =
    "/EDGE_NODE_NUM_ASSOCIATED_PARTITIONS_ARRAY";
const std::string EDGE_NODE_ASSOCIATED_PARTITIONS_DATASET =
    "/EDGE_NODE_ASSOCIATED_PARTITIONS_ARRAY";
const std::string DOSE_DATA_ARRAY = "/DOSE_ARRAY";
const std::string HIT_DATA_ARRAY = "/HIT_ARRAY";
const std::string DOSE_STATE_ARRAY = "/DOSE_STATE_ARRAY";
const std::string BEAM_DATASET = "/BEAM_DATA_ARRAY";

const std::string DOSE_RESULT_FOLDER_DATASET = "/DOSE_RESULT_FOLDER_NAME";
const std::string DOSE_NUM_CHUNKS_DATASET = "/DOSE_NUM_CHUNKS";
const std::string DOSE_NUM_VESSELS_PER_CHUNK_DATASET =
    "/DOSE_NUM_VESSELS_PER_CHUNK";
const std::string FIELD_SIZE_REDUCTION_ATTRIBUTE = "/FIELD_SIZE_REDUCTTION";