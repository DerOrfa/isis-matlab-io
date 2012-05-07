#include <DataStorage/io_interface.h>
#include <mclmcr.h>
#include <mat.h>
#include <boost/tuple/tuple.hpp>

namespace isis
{
namespace image_io
{
namespace _internal
{

class LibMatHandle{
public:
	static std::map<unsigned short, mxClassID> matID;
	static std::map<mxClassID,unsigned short> isisID;
	LibMatHandle();
    ~LibMatHandle();
    bool initialized;
};

class MatImageBase:public boost::noncopyable,public data::_internal::NDimensional<4>{
protected:
	mxArray *voxel_array,*root;
public:
	static LibMatHandle& handle();
	uint8_t *get(isis::util::FixedVector< size_t, 4 > offset, size_t bpv);
	~MatImageBase();
	mxArray* propMap2Array(const util::PropertyMap &pMap);
};

class MatImageRead:public MatImageBase{
	bool parseStruct(mxArray* branch, isis::util::PropertyMap& propParent, boost::shared_ptr< isis::data::Chunk > &chunk, std::string entry );
	boost::shared_ptr< data::Chunk > makeChunk();
	std::list<data::Chunk> &chunk_list;
public:
	MatImageRead(mxArray* src, std::string entry,std::list<data::Chunk> &target_chunk_list);
    static void sanitise(data::Chunk &chunk,size_t acNum);
	static bool hasOrTell(const char* name, isis::LogLevel log, const util::PropertyMap& pmap);
	template<typename T> static bool transformOrTell(const char* sname,const char* dname,LogLevel log,util::PropertyMap &pmap){
		return hasOrTell(sname,log,pmap) ?
			pmap.transform<T>(pmap.find(sname),dname):
			false;
	}
	template<typename T> static bool transformOrReplace(const char* sname,const char* dname,LogLevel log,util::PropertyMap &pmap,T replacement)
	{
		const util::PropertyMap::KeyType found=pmap.find(sname);
		util::checkType<T>();
		if(found.empty()){
			LOG(Runtime,log)
				<< "No property " << util::MSubject(sname) << " found, using default "
				<< util::Value<T>(replacement) << " for " << util::MSubject(dname);
			return pmap.setPropertyAs<T>(dname,replacement).isEmpty();
		} else if( pmap.propertyValue(found).isEmpty()) {
			LOG(Runtime,log)
				<< "Property " << util::MSubject(found) << " is empty, using default "
				<< util::Value<T>(replacement);
			pmap.remove(found);
			return pmap.setPropertyAs<T>(dname,replacement).isEmpty();
		} else
			return pmap.transform<T>(found,dname);
	}

};

class MatImageWrite:public MatImageBase{
	mxArray* prop2Array(const util::PropertyValue &prop);
	std::pair<unsigned short,mxClassID> getTypeForMat(const data::Image &img);
	mxArray* propMap2Array(const util::PropertyMap &pMap);
public:
	MatImageWrite(data::Image dat);
	bool good();
	bool write(MATFile* file,std::string entry_name);
};

class MatChunk:public data::Chunk{
	friend class MatImageRead;
    MatChunk(data::ValueArrayReference ref, util::FixedVector<size_t, 4 > size): Chunk(ref, size[1], size[0], size[2], size[3]){}
};

class Chunk2Array:public data::ChunkOp{
	MatImageWrite &m_target;
	size_t bpv;
public:
	Chunk2Array(MatImageWrite &target);
    virtual bool operator()(data::Chunk& , util::vector4< size_t > posInImage);
};


}}};