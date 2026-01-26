/*

Copyright (c) 2010--2011, Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the author at <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef __NABO_PRIVATE_H
#define __NABO_PRIVATE_H

#include "nabo.h"

#include <cstdint>
using std::uint32_t;

// Unused macro, add support for your favorite compiler
#if defined(__GNUC__)
	#define _UNUSED __attribute__ ((unused))
#else
	#define _UNUSED
#endif

/*!	\file nabo_private.h
	\brief header for implementation
	\ingroup private
*/

namespace Nabo
{
	//! \defgroup private private implementation
	//@{
	
	//! Euclidean distance
	template<typename T, typename A, typename B>
	inline T dist2(const A& v0, const B& v1)
	{
		return (v0 - v1).squaredNorm();
	}
	
	//! KDTree, unbalanced, points in leaves, stack, implicit bounds, ANN_KD_SL_MIDPT, optimised implementation
	template<typename T, typename Heap, typename CloudType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
	struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt: public NearestNeighbourSearch<T, CloudType>
	{
		typedef typename NearestNeighbourSearch<T, CloudType>::Vector Vector;
		typedef typename NearestNeighbourSearch<T, CloudType>::Matrix Matrix;
		typedef typename NearestNeighbourSearch<T, CloudType>::Index Index;
		typedef typename NearestNeighbourSearch<T, CloudType>::IndexVector IndexVector;
		typedef typename NearestNeighbourSearch<T, CloudType>::IndexMatrix IndexMatrix;
		
		using NearestNeighbourSearch<T, CloudType>::dim;
		using NearestNeighbourSearch<T, CloudType>::cloud;
		using NearestNeighbourSearch<T, CloudType>::creationOptionFlags;
		using NearestNeighbourSearch<T, CloudType>::minBound;
		using NearestNeighbourSearch<T, CloudType>::maxBound;
		using NearestNeighbourSearch<T, CloudType>::checkSizesKnn;
		
	protected:
		//! indices of points during kd-tree construction
		typedef std::vector<Index> BuildPoints;
		//! iterator to indices of points during kd-tree construction
		typedef typename BuildPoints::iterator BuildPointsIt;
		//! const-iterator to indices of points during kd-tree construction
		typedef typename BuildPoints::const_iterator BuildPointsCstIt;
		
		//! size of bucket
		const unsigned bucketSize;
		
		//! number of bits required to store dimension index + number of dimensions
		const uint32_t dimBitCount;
		//! mask to access dim
		const uint32_t dimMask;
		
		//! create the compound index containing the dimension and the index of the child or the bucket size
		inline uint32_t createDimChildBucketSize(const uint32_t dim, const uint32_t childIndex) const
		{ return dim | (childIndex << dimBitCount); }
		//! get the dimension out of the compound index
		inline uint32_t getDim(const uint32_t dimChildBucketSize) const
		{ return dimChildBucketSize & dimMask; }
		//! get the child index or the bucket size out of the coumpount index
		inline uint32_t getChildBucketSize(const uint32_t dimChildBucketSize) const
		{ return dimChildBucketSize >> dimBitCount; }
		
		//! search node
		struct Node
		{
			uint32_t dimChildBucketSize; //!< cut dimension for split nodes (dimBitCount lsb), index of right node or number of bucket(rest). Note that left index is current+1
			union
			{
				T cutVal; //!< for split node, split value
				uint32_t bucketIndex; //!< for leaf node, pointer to bucket
			};
			
			//! construct a split node
			Node(const uint32_t dimChild, const T cutVal):
				dimChildBucketSize(dimChild), cutVal(cutVal) {}
			//! construct a leaf node
			Node(const uint32_t bucketSize, const uint32_t bucketIndex):
				dimChildBucketSize(bucketSize), bucketIndex(bucketIndex) {}
		};
		//! dense vector of search nodes, provides better memory performances than many small objects
		typedef std::vector<Node> Nodes;
		
		//! search nodes
		Nodes nodes;
		
		//! buckets
		BuildPoints buckets;
		
		//! return the bounds of points from [first..last[ on dimension dim
		std::pair<T,T> getBounds(const BuildPointsIt first, const BuildPointsIt last, const unsigned dim);

		//! construct nodes for points [first..last[ inside the hyperrectangle [minValues..maxValues]
		unsigned buildNodes(BuildPointsIt first, BuildPointsIt last);
		
		//! search one point, call recurseKnn with the correct template parameters
		/** \param query pointer to query coordinates 
		 *	\param indices indices of nearest neighbours, must be of size k x query.cols()
		 *	\param dists2 squared distances to nearest neighbours, must be of size k x query.cols() 
		 *	\param i index of point to search
		 * 	\param heap reference to heap
		 * 	\param off reference to array of offsets
		 *	\param maxError error factor (1 + epsilon) 
		 *	\param maxRadius2 square of maximum radius
		 *	\param allowSelfMatch whether to allow self match
		 *	\param collectStatistics whether to collect statistics
		 *	\param sortResults wether to sort results
		 */
		unsigned long onePointKnn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, int i, Heap& heap, std::vector<T>& off, const T maxError, const T maxRadius2, const bool allowSelfMatch, const bool collectStatistics, const bool sortResults) const;
		
		//! recursive search, strongly inspired by ANN and [Arya & Mount, Algorithms for fast vector quantization, 1993]
		/**	\param query pointer to query coordinates 
		 * 	\param n index of node to visit
		 * 	\param rd squared dist to this rect
		 * 	\param heap reference to heap
		 * 	\param off reference to array of offsets
		 * 	\param maxError error factor (1 + epsilon) 
		 *	\param maxRadius2 square of maximum radius
		 */
		template<bool allowSelfMatch, bool collectStatistics>
		unsigned long recurseKnn(const T* query, const unsigned n, T rd, Heap& heap, std::vector<T>& off, const T maxError, const T maxRadius2) const;
		
	public:
		//! constructor, calls NearestNeighbourSearch<T>(cloud)
		KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags, const Parameters& additionalParameters);
		virtual unsigned long knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Index k, const T epsilon, const unsigned optionFlags, const T maxRadius) const;
		virtual unsigned long knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Vector& maxRadii, const Index k = 1, const T epsilon = 0, const unsigned optionFlags = 0) const;
	};
	
	//@}
}

#endif // __NABO_H
