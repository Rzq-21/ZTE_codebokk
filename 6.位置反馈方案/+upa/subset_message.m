function out=subset_message(message,refs)
[ok,slots]=ismember(refs,message.reference_ap_ids);assert(all(ok));out=message;out.reference_ap_ids=refs(:).';out.integer_codes=message.integer_codes(:,slots);
end
