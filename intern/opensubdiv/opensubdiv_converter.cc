/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2015 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sergey Sharybin.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include <cstdio>
#include <vector>

#ifdef _MSC_VER
#  include "iso646.h"
#endif

#include <opensubdiv/far/topologyRefinerFactory.h>

#include "opensubdiv_converter_capi.h"
#include "opensubdiv_intern.h"

#include <stack>

#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
namespace {

inline void reverse_face_verts(int *face_verts, int num_verts)
{
	int last_vert = face_verts[num_verts - 1];
	for (int i = num_verts - 1; i > 0;  --i) {
		face_verts[i] = face_verts[i - 1];
	}
	face_verts[0] = last_vert;
}

}  /* namespace */
#endif /* OPENSUBDIV_ORIENT_TOPOLOGY */

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {
namespace Far {

namespace {

template <typename T>
inline int findInArray(T array, int value)
{
	return (int)(std::find(array.begin(), array.end(), value) - array.begin());
}

#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
inline int get_loop_winding(int vert0_of_face, int vert1_of_face)
{
	int delta_face = vert1_of_face - vert0_of_face;
	if (abs(delta_face) != 1) {
		if (delta_face > 0) {
			delta_face = -1;
		}
		else {
			delta_face = 1;
		}
	}
	return delta_face;
}

inline void reverse_face_loops(IndexArray face_verts, IndexArray face_edges)
{
	for (int i = 0; i < face_verts.size() / 2; ++i) {
		int j = face_verts.size() - i - 1;
		if (i != j) {
			std::swap(face_verts[i], face_verts[j]);
			std::swap(face_edges[i], face_edges[j]);
		}
	}
	reverse_face_verts(&face_verts[0], face_verts.size());
}

inline void check_oriented_vert_connectivity(const int num_vert_edges,
                                             const int num_vert_faces,
                                             const int *vert_edges,
                                             const int *vert_faces,
                                             const int *dst_vert_edges,
                                             const int *dst_vert_faces)
{
#  ifndef NDEBUG
	for (int i = 0; i < num_vert_faces; ++i) {
		bool found = false;
		for (int j = 0; j < num_vert_faces; ++j) {
			if (vert_faces[i] == dst_vert_faces[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			assert(!"vert-faces connectivity ruined");
		}
	}
	for (int i = 0; i < num_vert_edges; ++i) {
		bool found = false;
		for (int j = 0; j < num_vert_edges; ++j) {
			if (vert_edges[i] == dst_vert_edges[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			assert(!"vert-edges connectivity ruined");
		}
	}
#  else
	(void)num_vert_edges;
	(void)num_vert_faces;
	(void)vert_edges;
	(void)vert_faces;
	(void)dst_vert_edges;
	(void)dst_vert_faces;
#  endif
}
#endif

}  /* namespace */

template <>
inline bool TopologyRefinerFactory<OpenSubdiv_Converter>::resizeComponentTopology(
        TopologyRefiner& refiner,
        const OpenSubdiv_Converter& conv)
{
	/* Faces and face-verts */
	const int num_faces = conv.get_num_faces(&conv);
	setNumBaseFaces(refiner, num_faces);
	for (int face = 0; face < num_faces; ++face) {
		const int num_verts = conv.get_num_face_verts(&conv, face);
		setNumBaseFaceVertices(refiner, face, num_verts);
	}
	/* Edges and edge-faces. */
	const int num_edges = conv.get_num_edges(&conv);
	setNumBaseEdges(refiner, num_edges);
	for (int edge = 0; edge < num_edges; ++edge) {
		const int num_edge_faces = conv.get_num_edge_faces(&conv, edge);
		setNumBaseEdgeFaces(refiner, edge, num_edge_faces);
	}
	/* Vertices and vert-faces and vert-edges/ */
	const int num_verts = conv.get_num_verts(&conv);
	setNumBaseVertices(refiner, num_verts);
	for (int vert = 0; vert < num_verts; ++vert) {
		const int num_vert_edges = conv.get_num_vert_edges(&conv, vert),
		          num_vert_faces = conv.get_num_vert_faces(&conv, vert);
		setNumBaseVertexEdges(refiner, vert, num_vert_edges);
		setNumBaseVertexFaces(refiner, vert, num_vert_faces);
	}
	return true;
}

template <>
inline bool TopologyRefinerFactory<OpenSubdiv_Converter>::assignComponentTopology(
	  TopologyRefiner& refiner,
        const OpenSubdiv_Converter& conv)
{
	using Far::IndexArray;
	/* Face relations. */
	const int num_faces = conv.get_num_faces(&conv);
	for (int face = 0; face < num_faces; ++face) {
		IndexArray dst_face_verts = getBaseFaceVertices(refiner, face);
		conv.get_face_verts(&conv, face, &dst_face_verts[0]);
		IndexArray dst_face_edges = getBaseFaceEdges(refiner, face);
		conv.get_face_edges(&conv, face, &dst_face_edges[0]);
	}
	/* Edge relations. */
	const int num_edges = conv.get_num_edges(&conv);
	for (int edge = 0; edge < num_edges; ++edge) {
		/* Edge-vertices */
		IndexArray dst_edge_verts = getBaseEdgeVertices(refiner, edge);
		conv.get_edge_verts(&conv, edge, &dst_edge_verts[0]);
		/* Edge-faces */
		IndexArray dst_edge_faces = getBaseEdgeFaces(refiner, edge);
		conv.get_edge_faces(&conv, edge, &dst_edge_faces[0]);
	}
#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
	/* Make face normals consistent. */
	bool *face_used = new bool[num_faces];
	memset(face_used, 0, sizeof(bool) * num_faces);
	std::stack<int> traverse_stack;
	int face_start = 0, num_traversed_faces = 0;
	/* Traverse all islands. */
	while (num_traversed_faces != num_faces) {
		/* Find first face of any untraversed islands. */
		while (face_used[face_start]) {
			++face_start;
		}
		/* Add first face to the stack. */
		traverse_stack.push(face_start);
		face_used[face_start] = true;
		/* Go over whole connected component. */
		while (!traverse_stack.empty()) {
			int face = traverse_stack.top();
			traverse_stack.pop();
			IndexArray face_edges = getBaseFaceEdges(refiner, face);
			ConstIndexArray face_verts = getBaseFaceVertices(refiner, face);
			for (int edge_index = 0; edge_index < face_edges.size(); ++edge_index) {
				const int edge = face_edges[edge_index];
				ConstIndexArray edge_faces = getBaseEdgeFaces(refiner, edge);
				if (edge_faces.size() != 2) {
					/* Can't make consistent normals for non-manifolds. */
					continue;
				}
				ConstIndexArray edge_verts = getBaseEdgeVertices(refiner, edge);
				/* Get winding of the reference face. */
				int vert0_of_face = findInArray(face_verts, edge_verts[0]),
				    vert1_of_face = findInArray(face_verts, edge_verts[1]);
				int delta_face = get_loop_winding(vert0_of_face, vert1_of_face);
				for (int edge_face = 0; edge_face < edge_faces.size(); ++edge_face) {
					int other_face = edge_faces[edge_face];
					/* Never re-traverse faces, only move forward. */
					if (face_used[other_face]) {
						continue;
					}
					IndexArray other_face_verts = getBaseFaceVertices(refiner,
					                                                  other_face);
					int vert0_of_other_face = findInArray(other_face_verts,
					                                      edge_verts[0]),
					    vert1_of_other_face = findInArray(other_face_verts,
					                                      edge_verts[1]);
					int delta_other_face = get_loop_winding(vert0_of_other_face,
					                                        vert1_of_other_face);
					if (delta_face * delta_other_face > 0) {
						IndexArray other_face_verts = getBaseFaceVertices(refiner,
						                                                  other_face),
						           other_face_edges = getBaseFaceEdges(refiner,
						                                               other_face);
						reverse_face_loops(other_face_verts,
						                   other_face_edges);
					}
					traverse_stack.push(other_face);
					face_used[other_face] = true;
				}
			}
			++num_traversed_faces;
		}
	}
#endif  /* OPENSUBDIV_ORIENT_TOPOLOGY */
	/* Vertex relations */
	const int num_verts = conv.get_num_verts(&conv);
	for (int vert = 0; vert < num_verts; ++vert) {

		/* Vert-Faces */
		IndexArray dst_vert_faces = getBaseVertexFaces(refiner, vert);
		int num_vert_faces = conv.get_num_vert_faces(&conv, vert);
		int *vert_faces = new int[num_vert_faces];
		conv.get_vert_faces(&conv, vert, vert_faces);
		/* Vert-Edges */
		IndexArray dst_vert_edges = getBaseVertexEdges(refiner, vert);
		int num_vert_edges = conv.get_num_vert_edges(&conv, vert);
		int *vert_edges = new int[num_vert_edges];
		conv.get_vert_edges(&conv, vert, vert_edges);
#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
		/* ** Order vertex edges and faces in a CCW order. ** */
		memset(face_used, 0, sizeof(bool) * num_faces);
		/* Number of edges and faces added to the ordered array. */
		int edge_count_ordered = 0, face_count_ordered = 0;
		/* Add loose edges straight into the edges array. */
		bool has_fan_connections = false;
		for (int i = 0; i < num_vert_edges; ++i) {
			IndexArray edge_faces = getBaseEdgeFaces(refiner, vert_edges[i]);
			if (edge_faces.size() == 0) {
				dst_vert_edges[edge_count_ordered++] = vert_edges[i];
			}
			else if (edge_faces.size() > 2) {
				has_fan_connections = true;
			}
		}
		if (has_fan_connections) {
			/* OpenSubdiv currently doesn't give us clues how to handle
			 * fan face connections. and since handling such connections
			 * complicates the loop below we simply don't do special
			 * orientation for them.
			 */
			memcpy(&dst_vert_edges[0], vert_edges, sizeof(int) * num_vert_edges);
			memcpy(&dst_vert_faces[0], vert_faces, sizeof(int) * num_vert_faces);
			delete [] vert_edges;
			delete [] vert_faces;
			continue;
		}
		/* Perform at max numbder of vert-edges iteration and try to avoid
		 * deadlock here for malformed mesh.
		 */
		for (int global_iter = 0; global_iter < num_vert_edges; ++global_iter) {
			/* Numbr of edges and faces which are still to be ordered. */
			int num_vert_edges_remained = num_vert_edges - edge_count_ordered,
			    num_vert_faces_remained = num_vert_faces - face_count_ordered;
			if (num_vert_edges_remained == 0 && num_vert_faces_remained == 0) {
				/* All done, nothing to do anymore. */
				break;
			}
			/* Face, edge and face-vertex inndex to start traversal from. */
			int face_start = -1, edge_start = -1, face_vert_start = -1;
			if (num_vert_edges_remained == num_vert_faces_remained) {
				/* Vertex is eitehr complete manifold or is connected to seevral
				 * manifold islands (hourglass-like configuration), can pick up
				 * random edge unused and start from it.
				 */
				/* TODO(sergey): Start from previous edge from which traversal
				 * began at previous iteration.
				 */
				for (int i = 0; i < num_vert_edges; ++i) {
					face_start = vert_faces[i];
					if (!face_used[face_start]) {
						ConstIndexArray
						    face_verts = getBaseFaceVertices(refiner, face_start),
						    face_edges = getBaseFaceEdges(refiner, face_start);
						face_vert_start = findInArray(face_verts, vert);
						edge_start = face_edges[face_vert_start];
						break;
					}
				}
			}
			else {
				/* Special handle of non-manifold vertex. */
				for (int i = 0; i < num_vert_edges; ++i) {
					bool start_found = false;
					edge_start = vert_edges[i];
					IndexArray edge_faces = getBaseEdgeFaces(refiner, edge_start);
					if (edge_faces.size() == 1) {
						face_start = edge_faces[0];
						if (!face_used[face_start]) {
							ConstIndexArray
							    face_verts = getBaseFaceVertices(refiner, face_start),
							    face_edges = getBaseFaceEdges(refiner, face_start);
							face_vert_start = findInArray(face_verts, vert);
							if (edge_start == face_edges[face_vert_start]) {
								start_found = true;
								break;
							}
						}
					}
					if (start_found) {
						break;
					}
					/* Reset indices for sanity check below. */
					face_start = edge_start = face_vert_start =  -1;
				}
			}
			/* Sanity check. */
			assert(face_start != -1 &&
			       edge_start != -1 &&
			       face_vert_start != -1);
			/* Traverse faces starting from the current one. */
			int edge_first = edge_start;
			dst_vert_faces[face_count_ordered++] = face_start;
			dst_vert_edges[edge_count_ordered++] = edge_start;
			face_used[face_start] = true;
			while (edge_count_ordered < num_vert_edges) {
				IndexArray face_verts = getBaseFaceVertices(refiner, face_start);
				IndexArray face_edges = getBaseFaceEdges(refiner, face_start);
				int face_edge_start = face_vert_start;
				int face_edge_next = (face_edge_start > 0) ? (face_edge_start - 1) : (face_verts.size() - 1);
				Index edge_next = face_edges[face_edge_next];
				if (edge_next == edge_first) {
					/* Multiple manifolds found, stop for now and handle rest
					 * in the next iteration.
					 */
					break;
				}
				dst_vert_edges[edge_count_ordered++] = edge_next;
				if (face_count_ordered < num_vert_faces) {
					IndexArray edge_faces = getBaseEdgeFaces(refiner, edge_next);
					assert(edge_faces.size() != 0);
					if (edge_faces.size() == 1) {
						assert(edge_faces[0] == face_start);
						break;
					}
					else if (edge_faces.size() != 2) {
						break;
					}
					assert(edge_faces.size() == 2);
					face_start = edge_faces[(edge_faces[0] == face_start) ? 1 : 0];
					face_vert_start = findInArray(getBaseFaceEdges(refiner, face_start), edge_next);
					dst_vert_faces[face_count_ordered++] = face_start;
					face_used[face_start] = true;
				}
				edge_start = edge_next;
			}
		}
		/* Verify ordering doesn't ruin connectivity information. */
		assert(face_count_ordered == num_vert_faces);
		assert(edge_count_ordered == num_vert_edges);
		check_oriented_vert_connectivity(num_vert_edges,
		                                 num_vert_faces,
		                                 vert_edges,
		                                 vert_faces,
		                                 &dst_vert_edges[0],
		                                 &dst_vert_faces[0]);
		/* For the release builds we're failing mesh construction so instead
		 * of nasty bugs the unsupported mesh will simply disappear from the
		 * viewport.
		 */
		if (face_count_ordered != num_vert_faces ||
		    edge_count_ordered != num_vert_edges)
		{
			delete [] vert_edges;
			delete [] vert_faces;
			return false;
		}
#else  /* OPENSUBDIV_ORIENT_TOPOLOGY */
		memcpy(&dst_vert_edges[0], vert_edges, sizeof(int) * num_vert_edges);
		memcpy(&dst_vert_faces[0], vert_faces, sizeof(int) * num_vert_faces);
#endif  /* OPENSUBDIV_ORIENT_TOPOLOGY */
		delete [] vert_edges;
		delete [] vert_faces;
	}
#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
	delete [] face_used;
#endif
	populateBaseLocalIndices(refiner);
	return true;
};

template <>
inline bool TopologyRefinerFactory<OpenSubdiv_Converter>::assignComponentTags(
        TopologyRefiner& refiner,
        const OpenSubdiv_Converter& conv)
{
	typedef OpenSubdiv::Sdc::Crease Crease;

	int num_edges = conv.get_num_edges(&conv);
	for (int edge = 0; edge < num_edges; ++edge) {
		float sharpness;
		ConstIndexArray edge_faces = getBaseEdgeFaces(refiner, edge);
		if (edge_faces.size() == 2) {
			sharpness = conv.get_edge_sharpness(&conv, edge);
		}
		else {
			/* Non-manifold edges must be sharp. */
			sharpness = Crease::SHARPNESS_INFINITE;
		}
		setBaseEdgeSharpness(refiner, edge, sharpness);
	}

	/* OpenSubdiv expects non-manifold vertices to be sharp but at the
	 * time it handles correct cases when vertex is a corner of plane.
	 * Currently mark verts which are adjacent to a loose edge as sharp,
	 * but this decision needs some more investigation.
	 */
	int num_vert = conv.get_num_verts(&conv);
	for (int vert = 0; vert < num_vert; ++vert) {
		ConstIndexArray vert_edges = getBaseVertexEdges(refiner, vert);
		for (int edge_index = 0; edge_index < vert_edges.size(); ++edge_index) {
			int edge = vert_edges[edge_index];
			ConstIndexArray edge_faces = getBaseEdgeFaces(refiner, edge);
			if (edge_faces.size() == 0) {
				setBaseVertexSharpness(refiner, vert, Crease::SHARPNESS_INFINITE);
				break;
			}
		}
	}

	return true;
}

template <>
inline void TopologyRefinerFactory<OpenSubdiv_Converter>::reportInvalidTopology(
        TopologyError /*errCode*/,
        const char *msg,
        const OpenSubdiv_Converter& /*mesh*/)
{
	printf("OpenSubdiv Error: %s\n", msg);
}

}  /* namespace Far */
}  /* namespace OPENSUBDIV_VERSION */
}  /* namespace OpenSubdiv */

namespace {

OpenSubdiv::Sdc::SchemeType get_capi_scheme_type(OpenSubdiv_SchemeType type)
{
	switch(type) {
		case OSD_SCHEME_BILINEAR:
			return OpenSubdiv::Sdc::SCHEME_BILINEAR;
		case OSD_SCHEME_CATMARK:
			return OpenSubdiv::Sdc::SCHEME_CATMARK;
		case OSD_SCHEME_LOOP:
			return OpenSubdiv::Sdc::SCHEME_LOOP;
	}
	assert(!"Unknown sceme type passed via C-API");
	return OpenSubdiv::Sdc::SCHEME_CATMARK;
}

}  /* namespace */

struct OpenSubdiv_TopologyRefinerDescr *openSubdiv_createTopologyRefinerDescr(
        OpenSubdiv_Converter *converter)
{
	typedef OpenSubdiv::Sdc::Options Options;

	using OpenSubdiv::Far::TopologyRefinerFactory;
	OpenSubdiv::Sdc::SchemeType scheme_type =
	        get_capi_scheme_type(converter->get_type(converter));
	Options options;
	options.SetVtxBoundaryInterpolation(Options::VTX_BOUNDARY_EDGE_ONLY);
	options.SetCreasingMethod(Options::CREASE_UNIFORM);
	options.SetFVarLinearInterpolation(Options::FVAR_LINEAR_ALL);

	TopologyRefinerFactory<OpenSubdiv_Converter>::Options
	        topology_options(scheme_type, options);
#ifdef OPENSUBDIV_VALIDATE_TOPOLOGY
	topology_options.validateFullTopology = true;
#endif
	/* We don't use guarded allocation here so we can re-use the refiner
	 * for GL mesh creation directly.
	 */
	return (struct OpenSubdiv_TopologyRefinerDescr*)
	        TopologyRefinerFactory<OpenSubdiv_Converter>::Create(
	                *converter,
	                topology_options);
}

void openSubdiv_deleteTopologyRefinerDescr(
        OpenSubdiv_TopologyRefinerDescr *topology_refiner)
{
	delete (OpenSubdiv::Far::TopologyRefiner *)topology_refiner;
}

int openSubdiv_topologyRefinerGetSubdivLevel(
        const OpenSubdiv_TopologyRefinerDescr *topology_refiner)
{
	using OpenSubdiv::Far::TopologyRefiner;
	const TopologyRefiner *refiner = (const TopologyRefiner *)topology_refiner;
	return refiner->GetMaxLevel();
}

int openSubdiv_topologyRefinerGetNumVerts(
        const OpenSubdiv_TopologyRefinerDescr *topology_refiner)
{
	using OpenSubdiv::Far::TopologyLevel;
	using OpenSubdiv::Far::TopologyRefiner;
	const TopologyRefiner *refiner = (const TopologyRefiner *)topology_refiner;
	const TopologyLevel &base_level = refiner->GetLevel(0);
	return base_level.GetNumVertices();
}

int openSubdiv_topologyRefinerGetNumEdges(
        const OpenSubdiv_TopologyRefinerDescr *topology_refiner)
{
	using OpenSubdiv::Far::TopologyLevel;
	using OpenSubdiv::Far::TopologyRefiner;
	const TopologyRefiner *refiner = (const TopologyRefiner *)topology_refiner;
	const TopologyLevel &base_level = refiner->GetLevel(0);
	return base_level.GetNumEdges();
}

int openSubdiv_topologyRefinerGetNumFaces(
        const OpenSubdiv_TopologyRefinerDescr *topology_refiner)
{
	using OpenSubdiv::Far::TopologyLevel;
	using OpenSubdiv::Far::TopologyRefiner;
	const TopologyRefiner *refiner = (const TopologyRefiner *)topology_refiner;
	const TopologyLevel &base_level = refiner->GetLevel(0);
	return base_level.GetNumFaces();
}

int openSubdiv_topologyRefnerCompareConverter(
        const OpenSubdiv_TopologyRefinerDescr *topology_refiner,
        OpenSubdiv_Converter *converter)
{
	using OpenSubdiv::Far::ConstIndexArray;
	using OpenSubdiv::Far::TopologyRefiner;
	using OpenSubdiv::Far::TopologyLevel;
	const TopologyRefiner *refiner = (const TopologyRefiner *)topology_refiner;
	const TopologyLevel &base_level = refiner->GetLevel(0);
	const int num_verts = base_level.GetNumVertices();
	const int num_edges = base_level.GetNumEdges();
	const int num_faces = base_level.GetNumFaces();
	/* Quick preliminary check. */
	OpenSubdiv::Sdc::SchemeType scheme_type =
	        get_capi_scheme_type(converter->get_type(converter));
	if (scheme_type != refiner->GetSchemeType()) {
		return false;
	}
	if (converter->get_num_verts(converter) != num_verts ||
	    converter->get_num_edges(converter) != num_edges ||
	    converter->get_num_faces(converter) != num_faces)
	{
		return false;
	}
	/* Compare all edges. */
	for (int edge = 0; edge < num_edges; ++edge) {
		ConstIndexArray edge_verts = base_level.GetEdgeVertices(edge);
		int conv_edge_verts[2];
		converter->get_edge_verts(converter, edge, conv_edge_verts);
		if (conv_edge_verts[0] != edge_verts[0] ||
		    conv_edge_verts[1] != edge_verts[1])
		{
			return false;
		}
	}
	/* Compare all faces. */
	std::vector<int> conv_face_verts;
	for (int face = 0; face < num_faces; ++face) {
		ConstIndexArray face_verts = base_level.GetFaceVertices(face);
		if (face_verts.size() != converter->get_num_face_verts(converter,
		                                                       face))
		{
			return false;
		}
		conv_face_verts.resize(face_verts.size());
		converter->get_face_verts(converter, face, &conv_face_verts[0]);
		bool direct_match = true;
		for (int i = 0; i < face_verts.size(); ++i) {
			if (conv_face_verts[i] != face_verts[i]) {
				direct_match = false;
				break;
			}
		}
		if (!direct_match) {
			/* If face didn't match in direct direction we also test if it
			 * matches in reversed direction. This is because conversion might
			 * reverse loops to make normals consistent.
			 */
#ifdef OPENSUBDIV_ORIENT_TOPOLOGY
			reverse_face_verts(&conv_face_verts[0], conv_face_verts.size());
			for (int i = 0; i < face_verts.size(); ++i) {
				if (conv_face_verts[i] != face_verts[i]) {
					return false;
				}
			}
#else
			return false;
#endif
		}
	}
	/* Compare sharpness. */
	for (int edge = 0; edge < num_edges; ++edge) {
		ConstIndexArray edge_faces = base_level.GetEdgeFaces(edge);
		float sharpness = base_level.GetEdgeSharpness(edge);
		float conv_sharpness;
		if (edge_faces.size() == 2) {
			conv_sharpness = converter->get_edge_sharpness(converter, edge);
		}
		else {
			conv_sharpness = OpenSubdiv::Sdc::Crease::SHARPNESS_INFINITE;
		}
		if (sharpness != conv_sharpness) {
			return false;
		}
	}
	return true;
}
