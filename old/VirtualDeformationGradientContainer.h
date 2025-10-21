/**
 * @file VirtualDeformationGradientContainer.h
 * @brief Storage helpers for deformation gradients of multiple virtual fields.
 */
#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include "DeformationGradientField.h"

/**
 * @brief Aggregates deformation gradient histories for several virtual fields.
 */
class VirtualDeformationGradientCollection
{
public:
	struct Field
	{
		std::string id;
		DeformationGradientHistory history;
	};

	void Clear()
	{
		m_fields.clear();
		m_lookup.clear();
	}

	size_t Size() const { return m_fields.size(); }
	bool Empty() const { return m_fields.empty(); }

	Field& Add(const std::string& id)
	{
		if (!id.empty())
		{
			auto it = m_lookup.find(id);
			if (it != m_lookup.end())
			{
				m_fields[it->second].history.Clear();
				return m_fields[it->second];
			}
		}

		Field field;
		field.id = id;
		field.history.Clear();
		m_fields.push_back(std::move(field));
		const size_t index = m_fields.size() - 1;
		if (!id.empty()) m_lookup[id] = index;
		return m_fields.back();
	}

	Field& operator[](size_t index) { return m_fields.at(index); }
	const Field& operator[](size_t index) const { return m_fields.at(index); }

	Field& operator[](const std::string& id)
	{
		Field* f = Find(id);
		if (f == nullptr) throw std::out_of_range("virtual deformation gradient field not found: " + id);
		return *f;
	}

	const Field& operator[](const std::string& id) const
	{
		const Field* f = Find(id);
		if (f == nullptr) throw std::out_of_range("virtual deformation gradient field not found: " + id);
		return *f;
	}

	Field* Find(const std::string& id)
	{
		if (!id.empty())
		{
			auto it = m_lookup.find(id);
			if (it != m_lookup.end()) return &m_fields[it->second];
		}
		for (auto& field : m_fields)
		{
			if (field.id == id) return &field;
		}
		return nullptr;
	}

	const Field* Find(const std::string& id) const
	{
		if (!id.empty())
		{
			auto it = m_lookup.find(id);
			if (it != m_lookup.end()) return &m_fields[it->second];
		}
		for (const auto& field : m_fields)
		{
			if (field.id == id) return &field;
		}
		return nullptr;
	}

	std::vector<Field>& Data() { return m_fields; }
	const std::vector<Field>& Data() const { return m_fields; }

	using iterator = std::vector<Field>::iterator;
	using const_iterator = std::vector<Field>::const_iterator;

	iterator begin() { return m_fields.begin(); }
	iterator end() { return m_fields.end(); }
	const_iterator begin() const { return m_fields.begin(); }
	const_iterator end() const { return m_fields.end(); }
	const_iterator cbegin() const { return m_fields.cbegin(); }
	const_iterator cend() const { return m_fields.cend(); }

private:
	std::vector<Field> m_fields;
	std::unordered_map<std::string, size_t> m_lookup;
};
