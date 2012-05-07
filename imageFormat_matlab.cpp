#include <CoreUtils/singletons.hpp>
#include <CoreUtils/value.hpp>
#include "imageFormat_matlab.hpp"
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace isis
{
namespace image_io
{
namespace _internal{

template<typename C,typename T, typename A> std::basic_string<C,T,A> fixFieldname(std::basic_string<C,T,A> name){
	const boost::basic_regex<C> replace_list("[^[:word:]_]");
	const std::basic_string<C>  buff(name.c_str()); //regex only works with char_traits<C>
	std::basic_string<C>  ret=boost::regex_replace<boost::regex_traits<C>,C>(buff,replace_list,"_", boost::match_default | boost::format_all);
	if(ret[ret.length()-1]=='_')ret=ret.substr(0,ret.find_last_not_of("_"));
	return ret.c_str();
}

LibMatHandle::LibMatHandle():initialized(false){
		/* Initialize Matlab runtime */
		const char *oplist[]={"-nojvm"};
		if (!mclInitializeApplication((const char **)oplist,1)){
			LOG(Runtime,error) << "Could not initialize MATLAB interface";
		} else {
			LOG(Debug,info) << "MATLAB interface initialized ";
			initialized=true;
			matID[(uint16_t)data::ValueArray<  int8_t>::staticID]=  mxINT8_CLASS;
			matID[(uint16_t)data::ValueArray< uint8_t>::staticID]= mxUINT8_CLASS;
			matID[(uint16_t)data::ValueArray< int16_t>::staticID]= mxINT16_CLASS;
			matID[(uint16_t)data::ValueArray<uint16_t>::staticID]=mxUINT16_CLASS;
			matID[(uint16_t)data::ValueArray< int32_t>::staticID]= mxINT32_CLASS;
			matID[(uint16_t)data::ValueArray<uint32_t>::staticID]=mxUINT32_CLASS;
			matID[(uint16_t)data::ValueArray< int64_t>::staticID]= mxINT64_CLASS;
			matID[(uint16_t)data::ValueArray<uint64_t>::staticID]=mxUINT64_CLASS;
			matID[(uint16_t)data::ValueArray<   float>::staticID]=mxSINGLE_CLASS;
			matID[(uint16_t)data::ValueArray<  double>::staticID]=mxDOUBLE_CLASS;


			for(std::map<unsigned short, mxClassID>::iterator p=matID.begin();p!=matID.end();p++)
				isisID[p->second]=p->first;
		}
	}

LibMatHandle::~LibMatHandle(){
	LOG(Debug,info) << "Shuting down MATLAB interface";
// 	mclTerminateApplication(); //@todo blocks
}

std::map<unsigned short, mxClassID> LibMatHandle::matID;
std::map<mxClassID,unsigned short> LibMatHandle::isisID;

template<typename T> util::Value<T> array2val(const mxArray *dat){
	const T* data = (T*)mxGetData(dat);
	return util::Value<T>(*data);
}

template<typename S,typename D> util::Value<std::list<D> > array2list(const mxArray *dat,mwSize len){
	const S *data=(S*)mxGetData(dat);
	return util::Value<std::list<D> >(util::listToList<D>(data,data+len));
}

LibMatHandle& MatImageBase::handle(){return util::Singletons::get<_internal::LibMatHandle,0>();}

MatImageWrite::MatImageWrite(data::Image dat){
	// adapt the image
	std::pair< unsigned int, mxClassID > type=getTypeForMat(dat);
	dat.convertToType(type.first);//unify type of the image to a fitting one for matlab
	dat.spliceDownTo(data::sliceDim);//must be slices
	
	// create master structure
	const char *field_names[] = {"voxel", "info"};
	root=mxCreateStructMatrix(1,1, 2 , field_names);

	//create voxel buffer
	const unsigned short dims = dat.getRelevantDims();
	const mwSize size[4]={dat.getNrOfRows(),dat.getNrOfColumns(),dat.getNrOfSlices(),dat.getNrOfTimesteps()}; // row and column are swapped in matlab
	init(size);//set up the NDimensional

	// copy the voxels
	if((voxel_array=mxCreateNumericArray(dims,size,type.second,mxREAL))){
		Chunk2Array cp(*this);
		dat.foreachChunk(cp);
		mxSetField(root, 0, "voxel", voxel_array);
	} else	
		LOG(Runtime,error) << "could not create matlab array of shape " << size << " wont copy voxel data";

	// copy properties
 	mxSetField(root, 0, "info", propMap2Array(dat));
}

MatImageRead::MatImageRead(mxArray* src,std::string entry,std::list<data::Chunk> &target_chunk_list):chunk_list(target_chunk_list)
{
	root=src;
	isis::util::PropertyMap propRoot;
	boost::shared_ptr< isis::data::Chunk > chunk;
	parseStruct(root,propRoot,chunk,entry);
	if(chunk){
		LOG(Runtime,info) << "Found single image in " << util::MSubject(entry);
		chunk_list.push_back(*chunk);
	}
}

bool MatImageRead::parseStruct(mxArray* branch, isis::util::PropertyMap& propParent, boost::shared_ptr< isis::data::Chunk > &chunk, std::string entry )
{
	const mwSize *dim=mxGetDimensions(branch);
	boost::tuples::tuple< util::PropertyMap, util::PropertyValue, boost::shared_ptr< data::Chunk > > ret;
	const mxClassID class_id=mxGetClassID(branch);

	if(class_id==mxSTRUCT_CLASS){
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// its a sub-structure
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		int fields=mxGetNumberOfFields(branch);
		boost::shared_ptr< isis::data::Chunk > sub_chunk;
		util::PropertyMap sub_map;
		for(int i=0;i<fields;i++){
				if(! parseStruct(mxGetFieldByNumber(branch,0,i),sub_map,sub_chunk,mxGetFieldNameByNumber(branch,i))){
					LOG(Runtime,warning) << "Failed to parse sub-structure " << util::MSubject(mxGetFieldNameByNumber(branch,i));
				}
		}
		if(sub_chunk){ // ok seems like there was a chunk found on this level so use the propmap we made from the other entries for the chunk
			LOG(Debug,info) << "Found a chunk of size " << util::MSubject(sub_chunk->getSizeAsString())
				<< " and with " <<  sub_map.getFlatMap().size() << " properties in structure " << util::MSubject(entry);
			sub_chunk->join(sub_map);
			chunk_list.push_back(*sub_chunk);
		} else { // elsewise insert the propmap as branch of the parent map
			LOG(Debug,info) << "Found a property branch  with " << sub_map.getFlatMap().size() << " entries in structure " << util::MSubject(entry);
			propParent.branch(entry.c_str())=sub_map;
		}
		return (sub_chunk.get() || sub_map.getKeys().size());
	} else if(class_id==mxCHAR_CLASS){
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// its a text-field
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		assert(dim[0]==1); // @todo I guess this _could_ happen
		if(dim[1]){
			std::auto_ptr<char> buffer((char*)malloc(dim[1]+1));
			mxGetString(branch,buffer.get(),dim[1]+1);
			propParent.setPropertyAs<std::string>(entry.c_str(),std::string(buffer.get()));
			LOG(Debug,verbose_info) << "Parsed text entry " << util::MSubject(entry) << " as " << util::MSubject( buffer.get() );
			return true;
		} else return false;
	} else if(dim[0]>1){
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// its an image
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		LOG(Debug,info) << "Numeric field " << util::MSubject(entry) << " has more than one row, assuming its an image.";
		voxel_array=branch;
		chunk=makeChunk();
		return chunk.get()!=NULL;
	} else if(dim[1]==1){
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// its a numeric field
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		switch(class_id){
		case mxUINT8_CLASS: propParent.propertyValue(entry.c_str())=array2val<uint8_t >(branch);break;
		case mxUINT16_CLASS:propParent.propertyValue(entry.c_str())=array2val<uint16_t>(branch);break;
		case mxUINT32_CLASS:propParent.propertyValue(entry.c_str())=array2val<uint32_t>(branch);break;
		case mxINT8_CLASS:  propParent.propertyValue(entry.c_str())=array2val< int8_t >(branch);break;
		case mxINT16_CLASS: propParent.propertyValue(entry.c_str())=array2val< int16_t>(branch);break;
		case mxINT32_CLASS: propParent.propertyValue(entry.c_str())=array2val< int32_t>(branch);break;
		case mxSINGLE_CLASS:propParent.propertyValue(entry.c_str())=array2val<float >(branch);break;
		case mxDOUBLE_CLASS:propParent.propertyValue(entry.c_str())=array2val<double>(branch);break;
		default:
			LOG(Runtime,error) << util::MSubject(entry) << "is of Unsupported type " << mxGetClassName(branch);
			return false;
		}
		LOG_IF(propParent.hasProperty(entry.c_str()),Debug,verbose_info) << "Parsed numeric entry " << util::MSubject(entry) << " as " << propParent.propertyValue(entry.c_str());
		return true;
	} else {
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// its a numeric list
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		switch(class_id){
		case mxUINT8_CLASS: propParent.propertyValue(entry.c_str())=array2list<uint8_t ,int32_t>(branch,dim[1]);break;
		case mxUINT16_CLASS:propParent.propertyValue(entry.c_str())=array2list<uint16_t,int32_t>(branch,dim[1]);break;
		case mxUINT32_CLASS:propParent.propertyValue(entry.c_str())=array2list<uint32_t,int32_t>(branch,dim[1]);break;
		case mxINT8_CLASS:  propParent.propertyValue(entry.c_str())=array2list< int8_t ,int32_t>(branch,dim[1]);break;
		case mxINT16_CLASS: propParent.propertyValue(entry.c_str())=array2list< int16_t,int32_t>(branch,dim[1]);break;
		case mxINT32_CLASS: propParent.propertyValue(entry.c_str())=array2list< int32_t,int32_t>(branch,dim[1]);break;
		case mxSINGLE_CLASS:propParent.propertyValue(entry.c_str())=array2list< float,double>(branch,dim[1]);break;
		case mxDOUBLE_CLASS:propParent.propertyValue(entry.c_str())=array2list<double,double>(branch,dim[1]);break;
		default:
			LOG(Runtime,error) << entry << "is of Unsupported type " << mxGetClassName(branch);
			return false;
		}
		LOG_IF(propParent.hasProperty(entry.c_str()),Debug,verbose_info) << "Parsed numeric list " << util::MSubject(entry) << " as " << propParent.propertyValue(entry.c_str());
		return true;
	}
	return false;
}

void MatImageRead::sanitise(isis::data::Chunk& chunk,size_t acNum)
{
	transformOrReplace<uint16_t             >("sequenceNumber",     "sequenceNumber",     warning,chunk,acNum);
	transformOrReplace<uint32_t             >("acquisitionNumber",  "acquisitionNumber",  warning,chunk,acNum);

	transformOrReplace<util::fvector4       >("voxelSize",          "voxelSize",          warning,chunk,util::fvector4(1,1,1,0));
	transformOrReplace<util::fvector4       >("indexOrigin",        "indexOrigin",        warning,chunk,util::fvector4(0,0,0,0));
	transformOrReplace<util::fvector4       >("rowVec",             "rowVec",             warning,chunk,util::fvector4(1,0));
	transformOrReplace<util::fvector4       >("columnVec",          "columnVec",          warning,chunk,util::fvector4(0,1));
	
	transformOrTell<util::fvector4          >("sliceVec",           "sliceVec",           info, chunk);

	transformOrTell<boost::posix_time::ptime>("sequenceStart",      "sequenceStart",      info, chunk);
	transformOrTell<std::string             >("sequenceDescription","sequenceDescription",info, chunk);
	transformOrTell<std::string             >("subjectName",        "subjectName",        info, chunk);
	transformOrTell<uint16_t                >("subjectAge",         "subjectAge",         info, chunk);
	transformOrTell<boost::gregorian::date  >("subjectBirth",       "subjectBirth",       info, chunk);
	transformOrTell<uint16_t                >("subjectWeigth",      "subjectWeigth",      info, chunk);
	transformOrTell<util::fvector4          >("voxelGap",           "voxelGap",           info, chunk);
	transformOrTell<std::string             >("transmitCoil",       "transmitCoil",       info, chunk);
	transformOrTell<uint16_t                >("flipAngle",          "flipAngle",          info, chunk);
	transformOrTell<std::string             >("performingPhysician","performingPhysician",info, chunk);
	transformOrTell<uint16_t                >("numberOfAverages",   "numberOfAverages",   info, chunk);
	transformOrTell<util::fvector4          >("caPos",              "caPos",              info, chunk);
	transformOrTell<util::fvector4          >("cpPos",              "cpPos",              info, chunk);
	transformOrTell<util::fvector4          >("diffusionGradient",  "diffusionGradient",  info, chunk);
	transformOrTell<uint16_t                >("repetitionTime",     "repetitionTime",     info, chunk);
	transformOrTell<float                   >("echoTime",           "echoTime",           info, chunk);
	transformOrTell<float                   >("acquisitionTime",    "acquisitionTime",    info, chunk);
}

bool MatImageRead::hasOrTell(const char* name, LogLevel log,const util::PropertyMap& pmap)
{
	const util::PropertyMap::KeyType found=pmap.find(name);
	if(found.empty()){
		LOG(Runtime,log) << "Property " << name << " not found.";
		return false;
	} else
		return true;
}

std::pair<unsigned short,mxClassID> MatImageWrite::getTypeForMat(const data::Image &img)
{
	const std::map<unsigned short, mxClassID> &map=util::Singletons::get<_internal::LibMatHandle,0>().matID;
	const unsigned short type_id=img.getMajorTypeID();
	
	std::map<unsigned short, mxClassID>::const_iterator found=map.find(type_id);
	if(found!=map.end())
		return *found;
	else switch(type_id){
		case data::ValueArray<bool>::staticID:
			LOG(Runtime,info) << data::ValueArray<bool>::staticName() << " is not supported by matlab, falling back to " << data::ValueArray<uint8_t>::staticName();
			return std::make_pair(data::ValueArray<uint8_t>::staticID,mxUINT8_CLASS);
			break;
		default:
			LOG(Runtime,warning) << "Sorry the type " << util::getTypeMap()[type_id] << " is not supported by matlab, will try double";
			return std::make_pair(data::ValueArray<double>::staticID,mxDOUBLE_CLASS);
	}
}

mxArray* MatImageWrite::prop2Array(const util::PropertyValue &prop)
{
	mxArray *ret;
	if((*prop).isFloat() || (*prop).isInteger())
		ret = mxCreateDoubleScalar(prop.as<double>());
	else switch(prop.getTypeID()){
		case util::Value<util::ivector4>::staticID:
		case util::Value<util::fvector4>::staticID:
		case util::Value<util::dvector4>::staticID:
			ret=mxCreateDoubleMatrix(1, 4, mxREAL);
			prop.as<util::dvector4>().copyTo(static_cast<double*>(mxGetData(ret)));//copy elements as double into the array
			break;
		case util::Value<util::ilist>::staticID:
		case util::Value<util::dlist>::staticID:{
			const util::dlist buff=prop.as<util::dlist>();
			ret=mxCreateDoubleMatrix(1, buff.size(), mxREAL);
			std::copy(buff.begin(),buff.end(),static_cast<double*>(mxGetData(ret)));//copy elements as double into the array
			}break;
		default:
			LOG(Debug,info) << "Falling back to string for " << prop.toString(true);
		case util::Value<std::string>::staticID:
		case util::Value<util::istring>::staticID:
			ret = mxCreateString(prop.toString(true).c_str());
			break;
	}
	LOG_IF(ret==NULL,Runtime,error) << "Failed to create matlab field from " << prop;

// 	mxSetField(branch, idx, name.c_str(),mxCreateString(prop.toString().c_str()));
// 	mxSetField(root, index, "seqNum", mxCreateDoubleScalar(seqNum));
// 	mxSetField(root, index, "center",vector2Array(pair.first.geometry.get_center()));
// 	mxSetField(root, index, "rotation",rotMatrix2Array(pair.first.geometry.get_gradrotmatrix()));
	return ret;
}

mxArray* MatImageWrite::propMap2Array(const util::PropertyMap& pMap)
{
	const util::PropertyMap::KeyList keys=pMap.getKeys();
	util::PropertyMap::KeyList root;
	BOOST_FOREACH(util::PropertyMap::KeyList::const_reference val,keys){
		root.insert(val.substr(0,val.find_first_of('/')));
	}
	static const mwSize one[]={1,1};

	mxArray *ret=mxCreateStructArray(2, one, 0, NULL);

	BOOST_FOREACH(const util::PropertyMap::KeyType &ref,root){
		util::istring name=fixFieldname(ref);
		const int fnum = mxAddField(ret, name.c_str());
		
		if(fnum>=0){
			if(pMap.hasBranch(ref)){ // if its a subtree
				mxSetFieldByNumber(ret, 0, fnum,propMap2Array(pMap.branch(ref))); // recurse
			} else {
				mxSetFieldByNumber(ret, 0, fnum,prop2Array(pMap.propertyValue(ref)));
			}
		} else
			LOG(Runtime,error) << "Failed to insert field " << util::MSubject( ref ) << " into the mat file";
	}
	return ret;
}


bool MatImageWrite::good(){return voxel_array!=NULL && root!=NULL;};

uint8_t *MatImageBase::get(util::FixedVector< size_t, 4 > offset,size_t bpv){
	const size_t idx[]={offset[data::columnDim],offset[data::rowDim],offset[data::sliceDim],offset[data::timeDim]}; // row and column are swapped in matlab
	const size_t at=getLinearIndex(idx);
	return static_cast<uint8_t*>(mxGetData(voxel_array))+at*bpv;
}

bool MatImageWrite::write(MATFile* file,std::string entry_name){
	return matPutVariable(file,entry_name.c_str(),root)==0;
}

boost::shared_ptr<data::Chunk> MatImageRead::makeChunk()
{
	const mxClassID class_id=mxGetClassID(voxel_array);
	const unsigned short type_id=handle().isisID[class_id];
	const mwSize *dim=mxGetDimensions(voxel_array);
	util::FixedVector<size_t,4> size;
	data::ValueArrayReference ref;
	boost::shared_ptr<data::Chunk> ret;

	if(type_id){
		size.fill(1);size.copyFrom(dim,dim+mxGetNumberOfDimensions(voxel_array));
		init(size); // init image with the size of the matlab array - so rows and columns are swapped (columns are the fastest incrementing index)
		ref=data::ValueArrayBase::createByID(type_id,size.product());

		ret.reset(new MatChunk(ref,size));
		boost::shared_ptr<uint8_t> dst= boost::static_pointer_cast<uint8_t>(ret->asValueArrayBase().getRawAddress());
		uint8_t *dPtr=dst.get();
		size_t bpv=ret->getBytesPerVoxel();

		util::FixedVector<size_t,4> index;
		for(index[data::timeDim]=0;index[data::timeDim]<size[data::timeDim];index[data::timeDim]++)
			for(index[data::sliceDim]=0;index[data::sliceDim]<size[data::sliceDim];index[data::sliceDim]++)
				for(index[data::columnDim]=size[data::rowDim];index[data::columnDim]>0;){
					index[data::columnDim]--;
					for(index[data::rowDim]=size[data::columnDim];index[data::rowDim]>0;){
						index[data::rowDim]--;
						memcpy(dPtr, get(index,bpv), bpv);dPtr+=bpv;
					}
				}
	} else {
		FileFormat::throwGenericError(std::string("Type \"")+mxGetClassName(voxel_array)+" not supported");
	}
	return ret;
}

MatImageBase::~MatImageBase(){
	if(root)
		mxDestroyArray(root);
}

Chunk2Array::Chunk2Array(MatImageWrite &target):m_target(target){}

bool Chunk2Array::operator()(data::Chunk& data, util::vector4< size_t > posInImage){
	assert(data.getRelevantDims()==2); // has to be slice-wise
	assert(posInImage[data::rowDim]==0 && posInImage[data::columnDim]==0);
	const size_t bpv=data.getBytesPerVoxel();
	util::FixedVector<size_t,4> chSize=data.getSizeAsVector();
	boost::shared_ptr<uint8_t> source= boost::static_pointer_cast<uint8_t>(data.asValueArrayBase().getRawAddress());
	uint8_t *dest=m_target.get(posInImage,bpv);

	for(size_t x=0;x<chSize[data::rowDim];x++)
		for(size_t y=0;y<chSize[data::columnDim];y++){
					const size_t ch_idx[]={x,y,0,0};
					const size_t ch_at=data.getLinearIndex(ch_idx);
					memcpy(dest,source.get()+ch_at*bpv,bpv);dest+=bpv;
				}
	return true;
}


}

class ImageFormat_matlab: public FileFormat
{
protected:
    util::istring suffixes(io_modes /*modes = both*/)const{return "mat";}

public:
	std::string getName()const {return "matlab data in/output";}

	int load( std::list<data::Chunk> &chunks, const std::string &filename, const util::istring &/*dialect */)  throw( std::runtime_error & ) {
		if(!util::Singletons::get<_internal::LibMatHandle,0>().initialized){
			LOG(Runtime,error) << "Matlab interface is not initialized, cannot read anything";
			return 0;
		}
		MATFile *file=matOpen(filename.c_str(),"r");
		if(file){
			const char *name;
			mxArray *matArray=matGetNextVariable(file, &name);
			_internal::MatImageRead img(matArray,name,chunks);
			size_t acNum=0;
			BOOST_FOREACH(data::Chunk &chRef,chunks){
				_internal::MatImageRead::sanitise(chRef,acNum++);
			}
		} else {
			throwGenericError(std::string("Could not open ")+filename+" for writing");
		}
		return 0;
	}
	void write( const data::Image &/*image*/, const std::string &/*filename*/, const util::istring &/*dialect */)throw( std::runtime_error & ) {} // not used
	void write(const std::list< data::Image >& images, const std::string& filename, const util::istring & /*dialect*/)throw( std::runtime_error & ) {
		if(!util::Singletons::get<_internal::LibMatHandle,0>().initialized){
			LOG(Runtime,error) << "Matlab interface is not initialized, cannot write anything";
			return;
		}
		if(images.size()==0){
			LOG(Runtime,info) << "Nothing to write to " << filename;
		}else if(images.size()>1)
			LOG(Runtime,notice) << "Storing " << images.size() << " images in " << filename;
		else
			LOG(Runtime,notice) << "Writing " << filename;

		MATFile *file=matOpen(filename.c_str(),"wz");
		if(file){
			BOOST_FOREACH(const data::Image &image,images){
				_internal::MatImageWrite mat(image);
				std::string image_name=std::string("S")+image.getPropertyAs<std::string>("sequenceNumber");
				if(image.hasProperty("sequenceDescription")) // matlab doesnt like that
					image_name+=std::string("_")+image.getPropertyAs<std::string>("sequenceDescription");

				if(!mat.good() || mat.write(file,_internal::fixFieldname(image_name))==false)
					throwGenericError(std::string("Failed to write to ") + filename +".");
			}
			matClose(file);
		} else {
			throwGenericError(std::string("Could not open ")+filename+" for writing");
		}
	}
	bool tainted()const {return false;}//internal plugins are not tainted
};
}
}
isis::image_io::FileFormat *factory()
{
	return new isis::image_io::ImageFormat_matlab();
}
