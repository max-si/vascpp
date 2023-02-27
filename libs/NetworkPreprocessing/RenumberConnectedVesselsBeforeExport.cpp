#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "VesselGenerator/ImportExportCommon.h"
//#include <vascular/core/ParallelSortTest.h>
//#include <vascular/core/Timsort.h>
#include "PreprocessorVessel.h"
#include "PartitionFunctions.h"

using std::cout;
using std::endl;
using std::map;
using std::string;
using std::unordered_map;
using std::vector;
;

struct ConnectedVesselInformation
{
    long long numConnectedVessels = 0;
    long long connectedVesselStartIndex = 0;
    long long numVesselsInPartition = 0;
    long long vesselStartIndex = 0;
    short partitionNumber = 0;
};

struct NetworkInformation
{
    vector<long long> counts;
    vector<long long> offsets;
    vector<int> rankNumVessels;
};

NetworkInformation GetNetworkInformation(
    const std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    NetworkInformation output;
    output.counts.resize(mpiSize * mpiSize);
    output.offsets.resize(mpiSize * mpiSize);
    output.rankNumVessels.resize(mpiSize);

    //Get local counts
    vector<long long> localCounts(mpiSize);
    //localCounts.resize(mpiSize);

    for (auto& vessel : vessels)
    {
        ++localCounts[vessel.getPartitionNumber()];
    }

    //Gather all counts
    //double time = -omp_get_wtime();
    MPI_Allgather(localCounts.data(), mpiSize, MPI_LONG_LONG_INT,
        output.counts.data(), mpiSize, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
    // if (!mpiRank) printf("\tAllgather Time: %1.6e\n",time+omp_get_wtime());
    localCounts.clear();

    //CalculateTotalCounts
    // time=-omp_get_wtime();
    for (int i = 0; i < mpiSize; ++i)
        for (int j = 0; j < mpiSize; ++j)
            output.rankNumVessels[i] += output.counts[mpiSize * i + j];
    //if (!mpiRank) printf("\tTotalCount Time: %1.6e\n",time+omp_get_wtime());
    //time=-omp_get_wtime();
    //calculate Offsets
    for (int i = 0; i < mpiSize - 1; ++i)
    {
        output.offsets[i + 1] += output.offsets[i];
        for (int j = 0; j < mpiSize; ++j)
            output.offsets[i + 1] += output.counts[mpiSize * j + i];
        //--output.offsets[i + 1];
    }
    for (int i = 1; i < mpiSize; ++i)
        for (int j = 0; j < mpiSize; ++j)
        {
            output.offsets[j + i * mpiSize] +=
                output.counts[mpiSize * (i - 1) + j] +
                output.offsets[j + (i - 1) * mpiSize];
        }
    //if (!mpiRank) printf("\tOffset Time: %1.6e\n",time+omp_get_wtime());

    return std::move(output);
}

vector<long long> ExtractVesselIndexes(
    const std::vector<PreprocessorVessel>& vessels)
{
    vector<long long> output;
    output.reserve(vessels.size());

    for (auto& vessel : vessels)
    {
        output.push_back(vessel.getInitialArrayPosition());
    }

    return std::move(output);
}

std::unordered_map<long long, long long> InitializeMap(
    const std::vector<PreprocessorVessel>& vessels)
{
    unordered_map<long long, long long> output;

    for (auto& vessel : vessels)
        for (auto& index : vessel.getConnectedVessels())
        {
            if (output.find(index) == output.end())
                output.emplace(index, -1l);
        }

    return std::move(output);
}

void CommunicateIndexesInRing(vector<long long>& sendBuffer,
    vector<long long>& recieveBuffer, const int numToSend,
    const int numToRecieve)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    vector<MPI_Request> requests(2);
    recieveBuffer.resize(numToRecieve);
    if (mpiRank == 0)
    {
        MPI_Irecv(recieveBuffer.data(), numToRecieve, MPI_LONG_LONG_INT,
            mpiSize - 1, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(sendBuffer.data(), numToSend, MPI_LONG_LONG_INT, mpiRank + 1,
            0, MPI_COMM_WORLD, &requests[1]);
    }
    else if (mpiRank == mpiSize - 1)
    {
        MPI_Irecv(recieveBuffer.data(), numToRecieve, MPI_LONG_LONG_INT,
            mpiRank - 1, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(sendBuffer.data(), numToSend, MPI_LONG_LONG_INT, 0, 0,
            MPI_COMM_WORLD, &requests[1]);
    }
    else
    {
        MPI_Irecv(recieveBuffer.data(), numToRecieve, MPI_LONG_LONG_INT,
            mpiRank - 1, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(sendBuffer.data(), numToSend, MPI_LONG_LONG_INT, mpiRank + 1,
            0, MPI_COMM_WORLD, &requests[1]);
    }

    MPI_Waitall(2, requests.data(), MPI_STATUSES_IGNORE);
}

void ApplyMapConnectedVessels(const unordered_map<long long, long long>& map,
    vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    for (auto& vessel : vessels)
    {
        for (auto& x : vessel.getConnectedVessels())
        {
            x = map.find(x)->second;
        }
    }
}

void AddIndexesToMap(const vector<long long>& vesselIndexes,
    const NetworkInformation& netInfo, const int partIndex,
    unordered_map<long long, long long>& map)
{
    int mpiSize;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    long long arrayIndex = 0;
    for (int i = 0; i < mpiSize; ++i)
        for (int j = 0; j < netInfo.counts[mpiSize * partIndex + i]; ++j)
        {
            auto mapItr = map.find(vesselIndexes[arrayIndex++]);
            if (mapItr != map.end())
            {
                mapItr->second = j + netInfo.offsets[mpiSize * partIndex + i];
            }
        }
}

unordered_map<long long, long long> CreateIndexMap(const std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //double time=-omp_get_wtime();
    NetworkInformation networkInfo = GetNetworkInformation(vessels);

    vector<long long> vesselIndexes = ExtractVesselIndexes(vessels);
    //if (!mpiRank) printf("\tExtract Vessel Time: %1.6e\n",time+=omp_get_wtime());
    //time=-omp_get_wtime();
    unordered_map<long long, long long> indexMap = InitializeMap(vessels);
    //if (!mpiRank) printf("\tinit Map Time: %1.6e\n",time+=omp_get_wtime());
    //time=-omp_get_wtime();
    AddIndexesToMap(vesselIndexes, networkInfo, mpiRank, indexMap);
    // if (!mpiRank) printf("\tAdd Index to Map Time: %1.6e\n",time+=omp_get_wtime());
    //Loop through other processes
    vector<long long> sendBuffer(vesselIndexes.size());
    //double time1=0;
    //time=0;
    int recieveIndex, sendIndex;
    for (int i = 1; i < mpiSize; i++)
    {
        recieveIndex = (mpiRank - i + mpiSize) % mpiSize;
        sendIndex = (mpiRank - (i - 1) + mpiSize) % mpiSize;

        //Communicate Data
        sendBuffer.assign(vesselIndexes.begin(), vesselIndexes.end());
        //        time1=-omp_get_wtime();
        CommunicateIndexesInRing(sendBuffer, vesselIndexes,
            networkInfo.rankNumVessels[sendIndex],
            networkInfo.rankNumVessels[recieveIndex]);
        //       time1+=omp_get_wtime();
        //add to map
        //       time=-omp_get_wtime();
        AddIndexesToMap(vesselIndexes, networkInfo, recieveIndex, indexMap);
        //      time+=omp_get_wtime();
    }

    // if (!mpiRank) printf("\tAdd Index to Map Time: %1.6e\n",time);
    // if (!mpiRank) printf("\tCommunicate Index Time: %1.6e\n",time1);

    return std::move(indexMap);
};

unordered_map<long long, long long> CreateIndexMapBcast(
    const std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    NetworkInformation networkInfo = GetNetworkInformation(vessels);

    /*if (mpiRank==0) cout << "Network Info" << endl;
    for (int i = 0; i < mpiSize; i++)
    {
        if (i == mpiRank)
        {
            cout << "Rank: " << mpiRank << '\n';
            cout << "Counts" <<'\n';
            for (int k =0;k<mpiSize;++k) {
                for (int j = 0; j < mpiSize; ++j)
                    cout << networkInfo.counts[k*mpiSize+j]<<',';
                cout<<'\n';
            }
            cout << "Offsets" <<'\n';
            for (int k =0;k<mpiSize;++k) {
                for (int j = 0; j < mpiSize; ++j)
                    cout << networkInfo.offsets[k*mpiSize+j]<<',';
                cout<<'\n';
            }
            cout <<endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (mpiRank == 0) cout << endl;*/

    vector<long long> localVesselIndexes = ExtractVesselIndexes(vessels);
    unordered_map<long long, long long> indexMap = InitializeMap(vessels);

    //Loop through other processes
    vector<long long> buffer(localVesselIndexes.size());
    int recieveIndex, sendIndex;
    for (int i = 0; i < mpiSize; i++)
    {
        recieveIndex = (mpiRank - i + mpiSize) % mpiSize;
        sendIndex = (mpiRank - (i - 1) + mpiSize) % mpiSize;

        if (i == mpiRank)
            buffer.assign(localVesselIndexes.begin(), localVesselIndexes.end());
        else
            buffer.resize(networkInfo.rankNumVessels[i]);

        MPI_Bcast(buffer.data(), networkInfo.rankNumVessels[i],
            MPI_LONG_LONG_INT, i, MPI_COMM_WORLD);

        //add to map
        AddIndexesToMap(buffer, networkInfo, i, indexMap);
    }

    return std::move(indexMap);
};

void RenumberConnectedVesselsBeforeExport(std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    // double time1 = -omp_get_wtime();
    unordered_map<long long, long long> indexMap = CreateIndexMap(vessels);
    // if (!mpiRank) printf("Create Index Time: %1.6e \n",time1+omp_get_wtime());
    /*if (mpiRank==0) cout << "Map Printing" << endl;
    for (int i = 0; i < mpiSize; i++)
    {
        if (i == mpiRank)
        {
            cout << "Rank: " << mpiRank << '\n';
            for (auto& x:indexMap)
                cout  <<x.first <<",\t " <<x.second<<'\n';
            cout <<endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (mpiRank == 0) cout << endl;*/

    // time1 = -omp_get_wtime();
    ApplyMapConnectedVessels(indexMap, vessels);
    // if (!mpiRank) printf("Apply Index Time: %1.6e \n",time1+omp_get_wtime());
}